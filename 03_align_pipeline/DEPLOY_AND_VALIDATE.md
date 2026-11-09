# Deploy the build and validate it on one yellowhead sample

Validation = running the new pipeline on ONE sample the old pipeline already aligned, and comparing the new BAMs with the old ones from Dropbox, record by record. Identical means the new pipeline changes no alignment
(the gsdf calibration was made on the old recipe's alignments). Nothing is uploaded: the old BAMs on Dropbox stay untouched; the local folder of the sample is overwritten, which is fine.

## 1. Laptop: put the build in the repository

Use the EXACT file name of the zip you downloaded last (`ls -l ~/Downloads/m03_increment1_*.zip` shows what is there; a glob in the unzip command breaks when several are present):

    cd ~/Downloads && rm -rf m03_increment1 && unzip -q m03_increment1_20261109-1418.zip
    ls m03_increment1/03_align_pipeline | head -20                        # must list DEPLOY_AND_VALIDATE.md
    cd ~/CichlidSRSequencing && git status --short                        # nothing should be listed
    git checkout reorg_spec && git pull --ff-only origin reorg_spec      # brings the panel's schema 2 and the updated SPEC.md and decisions.md
    cp -R ~/Downloads/m03_increment1/03_align_pipeline/. 03_align_pipeline/
    python3 tools/run_tests.py 03                                         # expect: Ran 158 tests ... OK  (needs pandas; about 30 s). An older build says 146 or 154.
    git add -A && git status --short                                      # 27 new files (M instead of A for files of an earlier build that was already pushed)
    git commit -m "Module 03: align_pipeline: panel identity check, --panel-expect, pinned environment, commit read without git, deployment guide"
    git push origin reorg_spec

## 2. Yellowhead: a second checkout (the existing one is not touched)

Git is needed only to clone and pull, and it lives in the `genomics` environment (the `align` environment has none and does not need it: the pipeline reads the checkout's commit from the .git folder by itself).

    conda activate genomics
    git clone -b reorg_spec "$(git -C ~/CichlidSRSequencing remote get-url origin)" ~/CSRS_reorg       # the first time; later: cd ~/CSRS_reorg && git pull --ff-only origin reorg_spec
    cd ~/CSRS_reorg && git log -1 --format='%h %ad %s' --date=short                                  # must be the commit you just pushed from the laptop
    conda activate align && source tools/env.sh                                                      # the pinned environment: minimap2 2.31, gatk4 4.6.2.0, pysam 0.24.1, samtools 1.24
    python -m align_pipeline run --help | grep -c panel-expect                                       # must print 1 (0 = an older build: pull again)
    grep -c write_sidecar 03_align_pipeline/align_pipeline/panel.py                                  # must print 0 (1 = the first build)
    which samtools gatk minimap2 rclone                                                              # all four must be found
    python -c "import pysam, pandas, numpy, yaml, gspread, gspread_dataframe, requests; print('imports ok')"
    python -m sd_genotypes panel check --modules all                                                 # four lines, one definition key each

The keys at the start of the pilot are `gsdf 5d71ee8d6c2b`, `te 2af781e90119`, `lg5inv 80bbea57e10b`, `bchrom 75367d0dfba8`; different keys mean a definition changed: stop and ask.
The `align` environment is built from `03_align_pipeline/env.yml` (`conda env create -f 03_align_pipeline/env.yml`). Its minimap2, gatk4 and pysam versions are those of `genomics` (checked 2026-11-09), which is
what the validation needs.

## 3. Choose the sample (read-only)

    python - <<'PYEOF'
    import os
    from cichlid_common.file_manager import FileManager
    fm = FileManager(); fm.setGenome('Mzebra_GT3_NCBI'); fm.readAlignmentDatabase()
    d = fm.alignment_dt.dropna(how='all')
    local = set(os.listdir(fm.localBamRefDir)) if os.path.isdir(fm.localBamRefDir) else set()      # the samples whose BAMs are on this machine (the 45 of the yellowhead batch)
    x = d[d.SampleID.isin(local)].copy()
    if x.empty:
        raise SystemExit('none of the %d sample folders in %s has a row in the AlignmentDatabase: is this the machine that aligned them?' % (len(local), fm.localBamRefDir))
    x['Runs'] = x.RunIDs.astype(str).str.count(',,') + 1
    x['BamGB'] = x.BamSize.astype(float) / 1e9
    print(x.sort_values('TotalReads')[['SampleID', 'Runs', 'TotalReads', 'Coverage', 'BamGB']].head(6).round(1).to_string(index=False))
    one = x[x.Runs == 1].sort_values('TotalReads')
    print('\nsuggested: ' + (one if len(one) else x.sort_values('TotalReads')).iloc[0].SampleID)
    PYEOF

The proposal is **SAMEA1904833**: in the old run's log its alignment took 13.7 minutes against 30 to 50 for most samples and 130 to 144 for the largest, which fits the smallest sample of the batch (5.2x, about 50 million
reads). That is an inference from the timings; the command above is the check. Use the sample it suggests.

## 4. Validate (in tmux; nothing is uploaded)

    S=SAMEA1904833                                                        # or the sample step 3 suggested
    python -m align_pipeline run Mzebra_GT3_NCBI -s $S -r --list         # -r: the sample already has a row; it must be listed
    mkdir -p ~/old_ref/$S
    rclone copy "ptm_dropbox:/COS/BioSci/BioSci-McGrath/Apps/CichlidSequencingData/Bamfiles/Mzebra_GT3_NCBI/$S/" ~/old_ref/$S/ --include "$S.all.bam" --include "$S.chimeric.bam" --progress
    tmux new -s validate
    python -m align_pipeline run Mzebra_GT3_NCBI -s $S -r --no-upload --keep-local --samples-in-flight 1 --panel-expect 5d71ee8d6c2b 2af781e90119 80bbea57e10b 75367d0dfba8 2>&1 | tee ~/validate_$S.log

If the panel is not ready the run stops at once, before anything is aligned, and says why; add `--no-panel` to validate the alignment alone. Expect about an hour (an estimate: only the alignment time is known).
Then, in the same terminal:

    B=$(python -c "from cichlid_common.file_manager import FileManager as F; f = F(); f.setGenome('Mzebra_GT3_NCBI'); print(f.localBamRefDir)")$S
    python ~/CSRS_reorg/03_align_pipeline/validate_against_old.py ~/old_ref/$S $B $S         # all and chimeric must say IDENTICAL
    python -m sd_genotypes.panel.identity --json $B/$S.panel.json --bam $B/$S.all.bam; echo "exit $?"     # 0 = the JSON belongs to this BAM
    ls -l $B $B/slices
    S=$S python - <<'PYEOF'
    import json, os
    from cichlid_common.file_manager import FileManager
    S = os.environ['S']
    fm = FileManager(); fm.setGenome('Mzebra_GT3_NCBI'); fm.readAlignmentDatabase()
    row = fm.alignment_dt[fm.alignment_dt.SampleID == S].iloc[0]
    qc = json.load(open(fm.localBamRefDir + S + '/' + S + '.qc.json'))
    for k in ('TotalReads', 'Coverage'):
        print(f'{k:12s} new {float(qc[k]):>18,.4f}   old row {float(row[k]):>18,.4f}')
    print('DuplicateFraction', qc['DuplicateFraction'], ' MappedFraction', qc['MappedFraction'], ' QCWarnings', repr(qc['QCWarnings']))
    PYEOF

## 5. What counts as a pass, and what to send me

Pass: the run ends with `1 of 1 samples complete`; `all` and `chimeric` are IDENTICAL; `TotalReads` and `Coverage` equal the old row's (TotalReads counted by `samtools idxstats` here and by GATK there: they must match
exactly); the panel JSON exists, `panel.identity` exits 0, and `QCWarnings` is empty. Send me `~/validate_$S.log` (its last 60 lines are enough), the output of the `validate_against_old.py` line and of the last snippet, and `ls -l` of the folder.

## If it stops

* "The panel is not ready": the message says which definition is missing, or that the keys differ from the pilot's. Fix that, or rerun with `--no-panel` to validate the alignment.
* "not found on PATH": the environment lacks that tool (step 2).
* "No samples to run": forgot `-r`, or a typo in the sample name.
* A failed sample: `FAILED at <stage>` in the log, with the tail of the tool's error output; the temp folder of the sample has the full logs.
* Ctrl-C stops every command it started; run the same command again and it continues where it stopped.
