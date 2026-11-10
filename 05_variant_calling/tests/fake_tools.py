"""Fake bcftools, samtools and rclone for the build tests. They follow the command lines `variant_calling.build` issues and write small TEXT files (named .vcf.gz), so that the
orchestration (staging, calling per contig, resume, renaming, the manifest, the safety stops) can be tested. They say nothing about the real tools, which were not installed where this
was written: the first real run is the check of bcftools' behaviour (--ignore-RG, -G -, -a, -f GQ).

Switches (environment variables): FAKE_LOG (append one JSON line per call), FAKE_NO_GROUPS / FAKE_NO_IGNORE_RG (hide an option from the usage text), FAKE_FAIL_CONTIG + FAKE_FAIL_FLAG
(mpileup of that contig fails while the flag file exists), FAKE_EXTRA_COLUMNS (extra sample columns), FAKE_RECORDS (variants per contig, default 3)."""
import json, os, shutil, sys


def log(tool, argv):
    p = os.environ.get('FAKE_LOG')
    if p:
        with open(p, 'a') as f:
            f.write(json.dumps(dict(tool=tool, argv=argv)) + '\n')


def die(msg, code=1):
    sys.stderr.write(msg + '\n')
    return code


def read_vcf(path):
    head, recs = [], []
    for l in open(path).read().splitlines():
        (head if l.startswith('#') else recs).append(l)
    return head, recs


def samples_of(head):
    cols = [l for l in head if l.startswith('#CHROM')][0].split('\t')
    return cols[9:]


def write_vcf(path, head, recs):
    with open(path, 'w') as f:
        f.write('\n'.join(head + recs) + '\n')


def opt(argv, flag, n=1):
    i = argv.index(flag)
    return argv[i + 1] if n == 1 else argv[i + 1:i + 1 + n]


def bcftools(argv):
    log('bcftools', argv)
    if not argv:
        return die('Usage: bcftools [--version] <command>', 1)
    cmd = argv[0]
    if cmd == '--version':
        print('bcftools 1.99\nUsing htslib 1.99')
        return 0
    if cmd == 'call' and len(argv) == 1:
        sys.stderr.write('Usage: bcftools call [options] <in.vcf.gz>\n   -m, --multiallelic-caller\n   -P, --prior <float>\n' + ('' if os.environ.get('FAKE_NO_GROUPS') else '   -G, --group-samples <file>|-\n'))
        return 1
    if cmd == 'mpileup' and len(argv) == 1:
        sys.stderr.write('Usage: bcftools mpileup [options] in1.bam\n   -b, --bam-list FILE\n' + ('' if os.environ.get('FAKE_NO_IGNORE_RG') else '       --ignore-RG\n'))
        return 1
    if cmd == 'mpileup':
        paths = [l.strip() for l in open(opt(argv, '-b')) if l.strip()]
        region = opt(argv, '-r')
        contig = region.split(':')[0]
        flag = os.environ.get('FAKE_FAIL_FLAG')
        if os.environ.get('FAKE_FAIL_CONTIG') == contig and flag and os.path.exists(flag):
            return die(f'simulated failure of mpileup on {contig}')
        paths = paths + [f'extra{i}' for i in range(int(os.environ.get('FAKE_EXTRA_COLUMNS', '0')))]
        out = ['##fileformat=VCFv4.2', '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(paths)]
        recs = [f'{contig}\t{100 * (i + 1)}\t.\tA\tG\t50\t.\t.\tGT:DP\t' + '\t'.join('0/1:5' for _ in paths) for i in range(int(os.environ.get('FAKE_RECORDS', '3')))]
        sys.stdout.write('\n'.join(out + recs) + '\n')
        return 0
    if cmd == 'call':
        if '-G' not in argv or opt(argv, '-G') != '-':
            return die('call: -G - was expected here')
        data = sys.stdin.read()
        with open(opt(argv, '-o'), 'w') as f:
            f.write(data)
        return 0
    if cmd == 'index':
        if '-n' in argv:
            head, recs = read_vcf(argv[-1])
            print(len(recs))
            return 0
        f = argv[-1]
        if os.path.exists(f + '.tbi') and '-f' not in argv:
            return die('index exists')
        open(f + '.tbi', 'w').write('tbi')
        return 0
    if cmd == 'query' and argv[1] == '-l':
        head, _ = read_vcf(argv[2])
        print('\n'.join(samples_of(head)))
        return 0
    if cmd == 'concat':
        files = [l.strip() for l in open(opt(argv, '-f')) if l.strip()]
        head, allrecs = None, []
        for f in files:
            h, r = read_vcf(f)
            head = head or h
            allrecs += r
        write_vcf(opt(argv, '-o'), head, allrecs)
        return 0
    if cmd == 'reheader':
        names = [l.strip() for l in open(opt(argv, '-s')) if l.strip()]
        head, recs = read_vcf(argv[-1])
        head = [('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(names)) if l.startswith('#CHROM') else l for l in head]
        write_vcf(opt(argv, '-o'), head, recs)
        return 0
    if cmd in ('norm', '+fill-tags', 'view'):
        src = argv[1] if cmd == '+fill-tags' else argv[-1]          # `bcftools +fill-tags IN -Oz -o OUT -- -t TAGS`: the input is the first argument after the plugin
        shutil.copy(src, opt(argv, '-o'))
        return 0
    return die(f'fake bcftools: {cmd} is not implemented')


def samtools(argv):
    log('samtools', argv)
    if argv[:1] == ['--version']:
        print('samtools 1.99')
        return 0
    if argv[:1] == ['quickcheck']:
        return 1 if b'CORRUPT' in open(argv[-1], 'rb').read() else 0
    return die(f'fake samtools: {argv[:1]} is not implemented')


def rclone(argv):
    log('rclone', argv)
    if argv[:1] == ['copy']:
        src, dst = argv[1], argv[2]
        inc = [argv[i + 1] for i, a in enumerate(argv) if a == '--include']
        if not os.path.isdir(src):
            return die(f'directory not found: {src}', 3)
        os.makedirs(dst, exist_ok=True)
        for n in os.listdir(src):
            if not inc or n in inc:
                shutil.copy(os.path.join(src, n), os.path.join(dst, n))
        return 0
    if argv[:1] == ['lsf']:
        if os.path.exists(argv[1]):
            print(os.path.basename(argv[1]))
        return 0
    return die(f'fake rclone: {argv[:1]} is not implemented')


def main(tool, argv):
    return {'bcftools': bcftools, 'samtools': samtools, 'rclone': rclone}[tool](argv)
