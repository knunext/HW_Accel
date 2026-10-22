#!/usr/bin/env python3
"""
mini_assembler.py — 교육용 de Bruijn 그래프 조립기 (파이썬 표준 라이브러리만 사용)

  리드 읽기 → k-mer 세기 → 거르기 → 유니티그 만들기 → 팁 정리 → 콘티그 → N50

사용법
  # 1) 참조 유전체에서 가짜 리드(FASTQ) 만들기
  python mini_assembler.py simulate phix.fa reads.fq --cov 30 --len 150 --err 0.002

  # 2) 리드만 보고 조립하기 (--ref 를 주면 결과를 참조와 비교)
  python mini_assembler.py assemble reads.fq -k 31 --out contigs.fa --ref phix.fa

설계 원칙
  - 이해하기 쉬운 코드를 우선합니다. 빠른 코드는 2-bit 인코딩(encode_2bit)과
    PIM·SmartSSD 장에서 다룹니다.
  - k는 홀수로 씁니다. 홀수 k-mer는 자기 자신의 역상보와 같을 수 없어서
    회문(palindrome) 처리를 하지 않아도 됩니다.
"""
import argparse
import collections
import gzip
import random
import resource
import sys
import time

# ---------------------------------------------------------------- 1. 서열 기초
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    """역상보: 반대쪽 가닥을 5'→3' 방향으로 읽은 서열."""
    return s.translate(COMP)[::-1]


def canonical(kmer):
    """정규형 k-mer: k-mer와 그 역상보 중 사전순으로 작은 쪽.
    어느 가닥에서 읽었든 같은 열쇠(key)가 되게 합니다."""
    rc = revcomp(kmer)
    return kmer if kmer <= rc else rc


def read_fasta(path):
    """FASTA 파일의 모든 레코드를 이어 붙여 하나의 문자열로 돌려줍니다."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        return "".join(line.strip() for line in f if not line.startswith(">")).upper()


def read_fastq(path):
    """FASTQ는 4줄이 한 리드입니다: @이름 / 서열 / + / 품질. 서열 줄만 꺼냅니다."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                yield line.strip().upper()


# ---------------------------------------------------------------- 2. 리드 흉내내기
def simulate_reads(genome, coverage=30, read_len=150, err=0.002, seed=1):
    """샷건 시퀀싱 흉내: 무작위 위치·무작위 가닥에서 read_len 길이를 잘라 오고,
    염기마다 err 확률로 다른 염기로 바꿉니다(치환 오류만).
    돌려주는 값: (서열, 품질 문자열) 생성기"""
    rng = random.Random(seed)
    n_reads = int(coverage * len(genome) / read_len)
    G = len(genome)
    for _ in range(n_reads):
        start = rng.randrange(0, G - read_len + 1)
        r = genome[start:start + read_len]
        if rng.random() < 0.5:                       # 절반은 반대 가닥
            r = revcomp(r)
        r = list(r)
        qual = []
        for i in range(read_len):
            q = 40 - int(12 * i / read_len) - rng.randrange(0, 4)   # 3' 끝으로 갈수록 품질 하락
            if rng.random() < err:
                r[i] = rng.choice([b for b in "ACGT" if b != r[i]])
                q = rng.randrange(2, 15)                        # 오류 자리는 대개 품질이 낮음
            qual.append(chr(33 + q))                            # Phred+33
        yield "".join(r), "".join(qual)


def write_fastq(path, reads):
    with open(path, "w") as f:
        for i, (seq, qual) in enumerate(reads):
            f.write(f"@sim_{i} length={len(seq)}\n{seq}\n+\n{qual}\n")


# ---------------------------------------------------------------- 3. k-mer 세기
def count_kmers(reads, k):
    """모든 리드의 모든 k-mer를 정규형으로 세어 Counter로 돌려줍니다.
    역상보는 리드마다 한 번만 만들어 두고, 위치 i의 k-mer와 짝이 되는
    역상보 조각을 잘라 비교합니다(k-mer마다 revcomp를 부르면 느립니다)."""
    counts = collections.Counter()
    for r in reads:
        L = len(r)
        if L < k:
            continue
        rc = revcomp(r)
        for i in range(L - k + 1):
            fw = r[i:i + k]
            bw = rc[L - k - i:L - i]          # fw의 역상보
            km = fw if fw <= bw else bw
            if "N" not in km:
                counts[km] += 1
    return counts


def kmer_spectrum(counts, max_count=200):
    """k-mer 스펙트럼: '몇 번 나온 k-mer가 몇 종류인가'의 히스토그램."""
    hist = collections.Counter(min(c, max_count) for c in counts.values())
    return [hist.get(c, 0) for c in range(max_count + 1)]


def auto_threshold(spectrum):
    """스펙트럼의 첫 골짜기(오류 봉우리와 진짜 봉우리 사이 최저점)를 문턱으로 고릅니다."""
    for c in range(2, len(spectrum) - 1):
        if spectrum[c] <= spectrum[c + 1]:
            return max(2, c)
    return 2


# ---------------------------------------------------------------- 4. 그래프와 유니티그
def successors(kmer, kset):
    """kmer 뒤에 이어질 수 있는 k-mer들 (마지막 k-1 글자 + 염기 하나)."""
    suffix = kmer[1:]
    return [suffix + b for b in "ACGT" if canonical(suffix + b) in kset]


def predecessors(kmer, kset):
    """kmer 앞에 올 수 있는 k-mer들 (염기 하나 + 처음 k-1 글자)."""
    prefix = kmer[:-1]
    return [b + prefix for b in "ACGT" if canonical(b + prefix) in kset]


def extend_right(start, kset, used):
    """갈림길이 나올 때까지 오른쪽으로 한 염기씩 걷습니다."""
    path = [start]
    cur = start
    while True:
        nxt = successors(cur, kset)
        if len(nxt) != 1:                      # 막다른 길(0) 또는 갈림길(2+)
            break
        n = nxt[0]
        if len(predecessors(n, kset)) != 1:    # 합류 지점이면 멈춤
            break
        cn = canonical(n)
        if cn in used:                         # 한 바퀴 돈 원형 경로
            break
        used.add(cn)
        path.append(n)
        cur = n
    return path


def build_unitigs(kset, k):
    """압축 de Bruijn 그래프: 갈림길 없는 경로를 하나의 서열(유니티그)로 묶습니다.
    양쪽 가닥을 모두 다루기 위해 왼쪽 확장은 역상보에서 오른쪽으로 걸어 구합니다."""
    used = set()
    unitigs = []
    for seed in kset:
        if seed in used:
            continue
        used.add(seed)
        right = extend_right(seed, kset, used)
        left = extend_right(revcomp(seed), kset, used)
        seq_left = "".join(revcomp(x)[0] for x in reversed(left[1:]))
        seq_right = seed + "".join(x[-1] for x in right[1:])
        unitigs.append(seq_left + seq_right)
    return unitigs


def remove_tips(kset, unitigs, k, max_len=None):
    """팁(tip) 정리: 한쪽 끝이 막다른 길인 짧은 유니티그는 대개 리드 끝의 오류에서
    생긴 곁가지입니다. 그 k-mer들을 지웁니다. 지운 개수를 돌려줍니다."""
    max_len = max_len or 2 * k
    removed = 0
    for u in unitigs:
        if len(u) >= max_len:
            continue
        head, tail = u[:k], u[-k:]
        dead_left = len(predecessors(head, kset)) == 0
        dead_right = len(successors(tail, kset)) == 0
        if dead_left != dead_right:           # 한쪽만 막혀 있으면 팁
            for i in range(len(u) - k + 1):
                kset.discard(canonical(u[i:i + k]))
            removed += 1
    return removed


# ---------------------------------------------------------------- 5. 평가
def n50(lengths):
    """길이를 내림차순으로 더해 전체의 절반을 넘는 순간의 콘티그 길이."""
    total = sum(lengths)
    acc = 0
    for L in sorted(lengths, reverse=True):
        acc += L
        if acc * 2 >= total:
            return L
    return 0


def evaluate(contigs, reference):
    """콘티그가 참조의 어느 가닥에 그대로(오류 없이) 들어 있는지 셉니다."""
    rc_ref = revcomp(reference)
    exact = sum(1 for c in contigs if c in reference or c in rc_ref)
    return exact


# ---------------------------------------------------------------- 6. 2-bit 인코딩 (하드웨어로 가는 다리)
ENC = {"A": 0, "C": 1, "G": 2, "T": 3}


def encode_2bit(kmer):
    """A=00 C=01 G=10 T=11. k<=31이면 64비트 정수 하나에 들어갑니다.
    PIM의 DPU와 SmartSSD의 FPGA는 문자열 대신 이 정수를 다룹니다."""
    x = 0
    for b in kmer:
        x = (x << 2) | ENC[b]
    return x


def mix64(x):
    """64비트 해시 섞기(splitmix64 마무리 단계). 어느 DPU가 맡을지 정할 때 씁니다."""
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9 & 0xFFFFFFFFFFFFFFFF
    x = (x ^ (x >> 27)) * 0x94D049BB133111EB & 0xFFFFFFFFFFFFFFFF
    return x ^ (x >> 31)


# ---------------------------------------------------------------- 7. 전체 파이프라인
def assemble(reads_path, k=31, threshold=None, min_contig=None, ref=None, out=None, log=print):
    t = {}
    t0 = time.perf_counter()
    counts = count_kmers(read_fastq(reads_path), k)
    t["count"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    spec = kmer_spectrum(counts)
    thr = threshold or auto_threshold(spec)
    kset = {km for km, c in counts.items() if c >= thr}
    n_distinct = len(counts)
    del counts
    t["filter"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    unitigs = build_unitigs(kset, k)
    t["unitig"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    tips = remove_tips(kset, unitigs, k)
    unitigs = build_unitigs(kset, k) if tips else unitigs
    t["tips"] = time.perf_counter() - t0

    min_contig = min_contig or 2 * k
    contigs = sorted((u for u in unitigs if len(u) >= min_contig), key=len, reverse=True)
    lens = [len(c) for c in contigs]
    res = dict(k=k, threshold=thr, distinct_kmers=n_distinct, solid_kmers=len(kset),
               tips_removed=tips, n_unitigs=len(unitigs), n_contigs=len(contigs),
               total_len=sum(lens), largest=lens[0] if lens else 0, n50=n50(lens),
               times=t, peak_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
               spectrum=spec[:80])
    if ref:
        reference = read_fasta(ref)
        res["ref_len"] = len(reference)
        res["exact_contigs"] = evaluate(contigs, reference)
    if out:
        with open(out, "w") as f:
            for i, c in enumerate(contigs):
                f.write(f">contig_{i} len={len(c)}\n")
                for j in range(0, len(c), 80):
                    f.write(c[j:j + 80] + "\n")
    log(f"k={k} 문턱={thr} 서로 다른 k-mer={n_distinct:,} 남긴 k-mer={len(kset):,}")
    log(f"콘티그 {len(contigs):,}개, 총 {sum(lens):,} bp, 최대 {res['largest']:,} bp, N50 {res['n50']:,} bp")
    log("단계별 시간(s): " + ", ".join(f"{a} {b:.2f}" for a, b in t.items()))
    return res


def main():
    ap = argparse.ArgumentParser(description="교육용 de Bruijn 그래프 조립기")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("simulate", help="참조 FASTA에서 FASTQ 리드 만들기")
    s.add_argument("ref"); s.add_argument("out")
    s.add_argument("--cov", type=float, default=30); s.add_argument("--len", type=int, default=150)
    s.add_argument("--err", type=float, default=0.002); s.add_argument("--seed", type=int, default=1)
    a = sub.add_parser("assemble", help="FASTQ 리드 조립")
    a.add_argument("reads"); a.add_argument("-k", type=int, default=31)
    a.add_argument("-t", "--threshold", type=int, default=None)
    a.add_argument("--min-contig", type=int, default=None)
    a.add_argument("--ref", default=None); a.add_argument("--out", default=None)
    args = ap.parse_args()
    if args.cmd == "simulate":
        g = read_fasta(args.ref)
        write_fastq(args.out, simulate_reads(g, args.cov, args.len, args.err, args.seed))
    else:
        if args.k % 2 == 0:
            sys.exit("k는 홀수로 주세요 (예: 21, 31, 55)")
        assemble(args.reads, args.k, args.threshold, args.min_contig, args.ref, args.out)


if __name__ == "__main__":
    main()
