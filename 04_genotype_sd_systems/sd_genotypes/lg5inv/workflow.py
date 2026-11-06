"""One pool: bchrom's RunProcessor (download -> uBAM -> align -> markdup -> slice -> depth -> finish) with a slice stage that also keeps the unmapped reads."""
import os
from sd_genotypes.bchrom import workflow as W


class Lg5RunProcessor(W.RunProcessor):
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.unmapped_bam = os.path.join(os.path.dirname(self.p.slice_bam), self.id + '.unmapped.bam')

    def _want_unmapped(self):
        return bool((self.cfg.get('slices') or {}).get('unmapped'))

    def finished(self):
        if not super().finished():
            return False
        return not self._want_unmapped() or os.path.exists(self.unmapped_bam)

    def _slice(self):
        """The region slice (indexed) of bchrom, then every unmapped read (samtools flag 4, which includes unmapped mates placed at their mapped partner) in a second BAM. Both are
        written under .partial names and renamed; a re-run of the stage redoes both."""
        super()._slice()
        if not self._want_unmapped():
            return
        os.makedirs(os.path.dirname(self.unmapped_bam), exist_ok=True)
        part = self.unmapped_bam + '.partial.bam'
        if os.path.exists(part):
            os.remove(part)
        threads = str(max(1, int(self.cfg['resources'].get('mosdepth_threads', 4))))
        self.ex.run(['samtools', 'view', '-@', threads, '-b', '-f', '4', '-o', part, self.p.bam], 'unmapped slice')
        self.ex.run(['samtools', 'quickcheck', '-v', part], 'quickcheck unmapped slice')
        os.replace(part, self.unmapped_bam)
