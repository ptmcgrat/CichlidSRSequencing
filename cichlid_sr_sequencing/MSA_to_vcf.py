"""
MSA-based unique variant caller.

Given a multiple sequence alignment and the record id of a target sample,
identify SNVs, insertions, and deletions that appear only in that sample
and emit them with VCF-style REF/ALT alleles (anchor base on indels).

All reference-side bases (anchor, SNV REF, deleted bases) come from the
FASTA. The MSA's reference row is used only to track coordinates --
never to source REF bases -- because in these alignments N in the MSA
reference indicates that the multi-sample alignment was uncertain at
that column, not that the assembly has an N there. The target's MSA
row remains the source for the ALT side: inserted bases for INS and
the alt nucleotide for SNV.
"""

from Bio import AlignIO
from Bio import SeqIO

import pyfaidx
from helper_modules.file_manager import FileManager as FM

class MSAUniqueVariantCaller:
    """Identify variants unique to one sample in a multiple sequence alignment.

    A variant is "unique" when every non-target sequence agrees and the
    target differs. Three kinds of variants are reported:

    * SNV - target and all others have ACGT bases, target differs.
    * INS - all others have gaps, target has bases.
    * DEL - target has a gap, all others agree on the same ACGT base.

    Detection deliberately rejects columns where any non-target sample
    has an ambiguous (N) base, because in this caller's intended use
    case those N's flag alignment uncertainty and would otherwise
    produce spurious calls (e.g. REF=A, ALT=A SNVs from FASTA REF
    matching the target while the MSA reference shows N).

    To avoid fragmenting long indels when the alignment has minor noise
    inside the indel region, the state machine uses lenient continuation:

    * An INS run starts only when every non-target sequence has a gap
      AND the target has a base. Once open it stays open until the
      reference sample advances (its column has a base), so internal
      noise no longer splits long insertions.
    * A DEL run starts when every non-target sequence agrees on an ACGT
      base AND the target has a gap. It stays open until the target
      advances (target column has a base).

    Parameters
    ----------
    msa_file : str
        Path to the MSA file.
    reference_fasta : str
        Path to the reference FASTA (indexed with pyfaidx). This is the
        authoritative source for every REF-side base.
    sample_name : str
        MSA record id of the sample to call unique variants for.
    msa_format : str, optional
        Format passed to ``Bio.AlignIO.read`` (default ``"clustal"``).
    ref_sample_name : str, optional
        Record id of the sequence in the MSA used to track reference
        coordinates. If None, the first sequence that is not
        ``sample_name`` is used.
    ref_seq_name : str, optional
        Sequence name in the FASTA file. Defaults to ``ref_sample_name``.
    output_chromosome : str, optional
        Chromosome label emitted in the output. Defaults to ``ref_seq_name``.
    region_start : int, optional
        1-based coordinate in the reference FASTA of the first base of
        the reference sample's alignment (default 1). The original
        script used a 0-based start; to reproduce its coordinates pass
        ``MZ_start + 1``.
    comparison_samples : list[str], optional
        Restrict the "other" samples to this subset of record ids. If
        None (default), every sequence except ``sample_name`` is used.
    qual : int, optional
        Value placed in the QUAL column of the output (default 100).
    """

    HEADER = ["Chromosome", "Position", "Name", "Reference", "Alt",
              "Q", "Filter", "Info"]

    def __init__(self, msa_file, reference_fasta, sample_name,
                 msa_format="clustal", ref_sample_name=None,
                 ref_seq_name=None, output_chromosome=None,
                 region_start=1, comparison_samples=None, qual=100):
        self.msa_file = msa_file
        self.alignment = AlignIO.read(msa_file, msa_format)
        self.fasta = pyfaidx.Fasta(reference_fasta)
        self.sample_name = sample_name
        self.region_start = int(region_start)
        self.qual = qual

        ids = [rec.id for rec in self.alignment]
        self._index = {rid: i for i, rid in enumerate(ids)}

        if sample_name not in self._index:
            raise ValueError(
                f"Sample '{sample_name}' not found in MSA. "
                f"Available records: {ids}"
            )
        self.sample_idx = self._index[sample_name]

        if ref_sample_name is None:
            ref_sample_name = next(r for r in ids if r != sample_name)
        if ref_sample_name not in self._index:
            raise ValueError(
                f"Reference sample '{ref_sample_name}' not found in MSA. "
                f"Available records: {ids}"
            )
        self.ref_sample_name = ref_sample_name
        self.ref_idx = self._index[ref_sample_name]

        self.ref_seq_name = ref_seq_name or ref_sample_name
        self.output_chromosome = output_chromosome or self.ref_seq_name

        if comparison_samples is None:
            self.other_indices = [i for i in range(len(self.alignment))
                                  if i != self.sample_idx]
        else:
            missing = [s for s in comparison_samples if s not in self._index]
            if missing:
                raise ValueError(
                    f"Comparison samples not in MSA: {missing}"
                )
            self.other_indices = [self._index[s] for s in comparison_samples]

        if not self.other_indices:
            raise ValueError("Need at least one non-target sample to compare.")

    # ----------------------------------------------------------------- helpers
    def _fasta_slice(self, start_ref_pos, end_ref_pos):
        """FASTA bases for an inclusive 1-based reference-position range.

        ``start_ref_pos`` and ``end_ref_pos`` count reference bases from
        the start of the alignment region (1-based). Returns ``"N"`` for
        any positions that fall outside the available FASTA range.
        """
        if end_ref_pos < start_ref_pos:
            return ""
        start_fasta = self.region_start + start_ref_pos - 1
        end_fasta = self.region_start + end_ref_pos - 1
        length = end_ref_pos - start_ref_pos + 1
        if start_fasta < 1:
            return "N" * length
        try:
            seq = str(self.fasta[self.ref_seq_name][start_fasta - 1:end_fasta]).upper()
        except (KeyError, IndexError):
            return "N" * length
        # Pad if the FASTA was shorter than expected
        if len(seq) < length:
            seq = seq + "N" * (length - len(seq))
        return seq

    def _fasta_base(self, ref_pos):
        """Single FASTA base at the given 1-based reference position."""
        return self._fasta_slice(ref_pos, ref_pos) or "N"

    def validate_reference(self):
        """Verify the reference sample's MSA bases match the FASTA.

        Returns a list of ``(alignment_column, fasta_base, msa_base)``
        tuples for every non-gap MSA position that disagrees with the
        FASTA. Positions where the MSA shows ``N`` will appear as
        mismatches; that is expected when the MSA uses N to indicate
        alignment uncertainty.
        """
        mismatches = []
        ref_pos = 0
        aln = self.alignment
        for i in range(aln.get_alignment_length()):
            base = aln[self.ref_idx, i].upper()
            if base == "-":
                continue
            ref_pos += 1
            expected = self._fasta_base(ref_pos)
            if base != expected:
                mismatches.append((i, expected, base))
        return mismatches

    # --------------------------------------------------------- variant calling
    def call_variants(self):
        """Yield variant dictionaries.

        Yields
        ------
        dict
            Keys ``chromosome, position, name, ref, alt, qual, filter, info``.
        """
        aln = self.alignment
        s = self.sample_idx
        r = self.ref_idx
        others = self.other_indices

        ref_pos = 0           # 1-based count of reference-sample bases consumed
        state = None          # None | "INS" | "DEL"
        pending = None        # dict describing the currently-open indel

        for i in range(aln.get_alignment_length()):
            col = aln[:, i]
            tgt = col[s].upper()
            ref = col[r].upper()
            other_nucs = [col[j].upper() for j in others]

            if ref != "-":
                ref_pos += 1

            others_set = set(other_nucs)
            others_uniform = len(others_set) == 1
            others_consensus = other_nucs[0] if others_uniform else None
            others_all_gap = others_uniform and others_consensus == "-"
            others_all_acgt = (others_uniform
                              and others_consensus in "ACGT")

            is_ins_start = (tgt != "-") and others_all_gap
            is_del_start = (tgt == "-") and others_all_acgt
            is_snv = (others_all_acgt and tgt in "ACGT"
                      and tgt != others_consensus)

            # ---- Step 1: close any pending indel that has hit its boundary
            if state == "INS" and ref != "-":
                yield self._build_indel(pending)
                state = None
                pending = None
            elif state == "DEL" and tgt != "-":
                yield self._build_indel(pending)
                state = None
                pending = None

            # ---- Step 2: extend the pending indel, or start a new variant
            if state == "INS":
                # Inside an insertion region (ref still gap). Append any
                # target base for the ALT allele.
                if tgt != "-":
                    pending["alt_seq"] += tgt
                pending["col_end"] = i
            elif state == "DEL":
                # Inside a deletion region (target still gap). Track the
                # range of consumed reference positions; the actual REF
                # bases are pulled from the FASTA at emit time.
                if ref != "-":
                    pending["last_deleted_ref_pos"] = ref_pos
                pending["col_end"] = i
            else:
                if is_ins_start:
                    state = "INS"
                    pending = {
                        "kind": "INS",
                        "col_start": i,
                        "col_end": i,
                        "anchor_ref_pos": ref_pos,   # last consumed ref base
                        "alt_seq": tgt,              # from target sample
                    }
                elif is_del_start:
                    state = "DEL"
                    pending = {
                        "kind": "DEL",
                        "col_start": i,
                        "col_end": i,
                        # anchor is the ref base BEFORE the deleted one;
                        # ref_pos was just incremented for this column.
                        "anchor_ref_pos": ref_pos - 1,
                        "first_deleted_ref_pos": ref_pos,
                        "last_deleted_ref_pos": ref_pos,
                    }
                elif is_snv:
                    yield self._make_snv(i, ref_pos, tgt)

        # End of alignment: flush a trailing indel if one is still open.
        if state is not None:
            yield self._build_indel(pending)

    # ----------------------------------------------------- record constructors
    def _make_snv(self, col, ref_pos, alt_nuc):
        ref_nuc = self._fasta_base(ref_pos)
        return {
            "chromosome": self.output_chromosome,
            "position": self.region_start + ref_pos - 1,
            "name": ".",
            "ref": ref_nuc,
            "alt": alt_nuc,
            "qual": self.qual,
            "filter": ".",
            "info": f"TYPE=SNV;ALN_COL={col + 1}",
        }

    def _build_indel(self, p):
        anchor_pos = self.region_start + p["anchor_ref_pos"] - 1
        anchor_base = self._fasta_base(p["anchor_ref_pos"])

        if p["kind"] == "INS":
            alt_seq = p["alt_seq"]
            ref_allele = anchor_base
            alt_allele = anchor_base + alt_seq
            info = f"TYPE=INS;ALN_COL={p['col_start'] + 1};LEN={len(alt_seq)}"
        else:
            deleted = self._fasta_slice(
                p["first_deleted_ref_pos"], p["last_deleted_ref_pos"]
            )
            ref_allele = anchor_base + deleted
            alt_allele = anchor_base
            info = f"TYPE=DEL;ALN_COL={p['col_start'] + 1};LEN={len(deleted)}"

        return {
            "chromosome": self.output_chromosome,
            "position": anchor_pos,
            "name": ".",
            "ref": ref_allele,
            "alt": alt_allele,
            "qual": self.qual,
            "filter": ".",
            "info": info,
        }

    # --------------------------------------------------------------- I/O glue
    def write_tsv(self, output_file):
        """Write variants to a tab-delimited file matching the original layout."""
        with open(output_file, "w") as f:
            f.write("\t".join(self.HEADER) + "\n")
            for v in self.call_variants():
                f.write("\t".join(str(x) for x in (
                    v["chromosome"], v["position"], v["name"],
                    v["ref"], v["alt"], v["qual"], v["filter"], v["info"],
                )) + "\n")

    def write_vcf(self, output_file):
        """Write variants to a minimal but valid VCF 4.2 file."""
        with open(output_file, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write(f"##source=MSAUniqueVariantCaller\n")
            f.write(f"##sample={self.sample_name}\n")
            f.write(f"##reference={self.ref_seq_name}\n")
            f.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="SNV, INS, or DEL">\n')
            f.write('##INFO=<ID=ALN_COL,Number=1,Type=Integer,Description="1-based alignment column where the variant starts">\n')
            f.write('##INFO=<ID=LEN,Number=1,Type=Integer,Description="Length of inserted or deleted sequence">\n')
            f.write("#" + "\t".join([
                "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            ]) + "\n")
            for v in self.call_variants():
                f.write("\t".join(str(x) for x in (
                    v["chromosome"], v["position"], v["name"],
                    v["ref"], v["alt"], v["qual"],
                    "PASS" if v["filter"] == "." else v["filter"],
                    v["info"],
                )) + "\n")


if __name__ == "__main__":
    
    fm_obj = FM(genome_version = 'Mzebra_GT3_NCBI')
    count = SeqIO.convert("MappedRegionsMSA.aln", "clustal", "MappedRegionsMSA.fasta", "fasta")

    caller = MSAUniqueVariantCaller(
        msa_file="MappedRegionsMSA.aln",
        reference_fasta=fm_obj.localGenomeFile,
        sample_name="Y_reg",         # MSA record id of the target sample
        ref_sample_name="MZ_reg",           # MSA record id of the reference sample
        ref_seq_name="NC_135176.1",         # actual sequence name in the FASTA
        output_chromosome="NC_135176.1",
        region_start=24777955,              # 1-based first base; original used MZ_start=24777954
    )
    caller.write_tsv("UniqueVariantsMSA.csv")
    caller.write_vcf("UniqueVariantsMSA.vcf")