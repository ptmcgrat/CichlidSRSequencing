"""Disk-space accounting for the sequencing pipeline: a sample is only ADMITTED (its uBAMs downloaded) when the disk can hold everything it will write.

Why: a download queue that simply keeps going fills the disk when compute or uploads are slower than the network. Every sample in flight writes much more than its uBAMs:
the uBAMs, the sorted BAM, the duplicate-marked BAM, MarkDuplicates scratch files and the split BAMs, about 3 times the uBAM size at the peak (the old worker's own rule of thumb).

The ledger reserves `reserve_factor * (uBAM bytes)` for every admitted sample. A new sample is admitted only if

    free disk now  -  (what the samples already in flight have still to write)  >=  its own reservation + min_free

where "still to write" = reservation - what the sample already occupies on disk (measured). When there is not enough room the caller waits (a finished and cleaned-up sample frees
space). The reservation of a sample shrinks to its measured size once its large writes are over (settle), and disappears when its local files are deleted or it fails (release)."""
import os
import shutil
import threading


def dir_size(path):
    """Bytes in all files under path (0 if it does not exist)."""
    total = 0
    for root, _, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return total


class DiskLedger:
    def __init__(self, path, reserve_factor=5.5, min_free_gb=200, usage_fn=None, size_fn=None):
        """usage_fn(path) -> object with .free (bytes); size_fn(sample) -> bytes this sample currently occupies (uBAMs + outputs + scratch)."""
        self.path, self.factor, self.min_free = path, float(reserve_factor), float(min_free_gb) * 1e9
        self.usage_fn = usage_fn or shutil.disk_usage
        self.size_fn = size_fn or (lambda sample: 0)
        self.reserved, self.lock = {}, threading.Lock()

    def free(self):
        return self.usage_fn(self.path).free

    def _unwritten(self):
        return sum(max(0.0, r - self.size_fn(s)) for s, r in self.reserved.items())

    def admit(self, sample, need_bytes):
        """(True, '') and the reservation is made, or (False, reason). Admitting the same sample twice is harmless."""
        with self.lock:
            if sample in self.reserved:
                return True, ''
            want = self.factor * float(need_bytes)
            free, pending = self.free(), self._unwritten()
            if free - pending >= want + self.min_free:
                self.reserved[sample] = want
                return True, ''
            return False, (f'free disk {free / 1e9:.0f} GB, of which {pending / 1e9:.0f} GB are still to be written by samples in flight; this sample needs {want / 1e9:.0f} GB '
                           f'plus a {self.min_free / 1e9:.0f} GB margin')

    def settle(self, sample):
        """The sample's large writes are over: reserve only what it occupies now."""
        with self.lock:
            if sample in self.reserved:
                self.reserved[sample] = float(self.size_fn(sample))

    def release(self, sample):
        with self.lock:
            self.reserved.pop(sample, None)

    def idle(self):
        """No sample holds a reservation: nothing in flight will free any space."""
        with self.lock:
            return not self.reserved

    def snapshot(self):
        with self.lock:
            return dict(free_gb=self.free() / 1e9, reserved_samples=len(self.reserved), unwritten_gb=self._unwritten() / 1e9)
