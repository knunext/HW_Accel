"""realdata.py — 진짜 리드에 가까운 흉내와 정리, 버블 터뜨리기 (15장)

흉내: 짧은 조각(어댑터가 섞임), 3' 끝으로 갈수록 오류가 늘어나는 품질, GC 쏠림, 두 벌(이배체) 유전체
정리: 품질 창으로 끝 자르기, 어댑터 자르기
그래프: 이형접합 SNP가 만드는 버블을 찾아 약한 쪽을 지우기
"""
import math
import random

from mini_assembler import canonical, predecessors, revcomp, successors

ADAPTER = "AGATCGGAAGAGC"          # 일루미나 어댑터의 앞부분 (리드가 조각 끝을 넘어가면 나타남)


# ---------------------------------------------------------------- 흉내
def make_diploid(genome, snp_rate=0.005, seed=1):
    """두 번째 벌: snp_rate 비율로 염기를 바꾼 서열. 두 벌이 섞인 시료를 흉내 냅니다."""
    rng = random.Random(seed)
    h2 = list(genome)
    snps = []
    for i in range(len(genome)):
        if rng.random() < snp_rate:
            h2[i] = rng.choice([b for b in "ACGT" if b != genome[i]])
            snps.append(i)
    return "".join(h2), snps


def gc_weight(seq, center=0.45, width=0.12):
    """GC 비율이 center에서 멀수록 덜 읽히는 쏠림 (PCR 증폭 효율을 흉내)."""
    gc = (seq.count("G") + seq.count("C")) / max(1, len(seq))
    return math.exp(-((gc - center) / width) ** 2)


def simulate_real_reads(haplotypes, coverage=30, read_len=150, frag_mean=300, frag_sd=80,
                        err_start=0.001, err_end=0.02, gc_bias=True, seed=1):
    """조각 길이가 리드보다 짧으면 리드 끝에 어댑터가 붙습니다. 오류율은 3' 끝으로 갈수록 올라가고,
    품질 글자는 그 오류율에 맞춰 적습니다. 돌려주는 값: (서열, 품질) 생성기."""
    rng = random.Random(seed)
    G = len(haplotypes[0])
    n = int(coverage * G / read_len)
    made = 0
    while made < n:
        h = haplotypes[rng.randrange(len(haplotypes))]
        flen = max(40, int(rng.gauss(frag_mean, frag_sd)))
        s = rng.randrange(0, G - flen)
        frag = h[s:s + flen]
        if gc_bias and rng.random() > gc_weight(frag):
            continue                                        # 이 조각은 덜 증폭돼 안 읽힘
        if rng.random() < 0.5:
            frag = revcomp(frag)
        r = (frag + ADAPTER + "A" * read_len)[:read_len]     # 조각이 짧으면 어댑터로 넘어감
        seq, qual = [], []
        for i, b in enumerate(r):
            e = err_start + (err_end - err_start) * (i / (read_len - 1)) ** 2
            if rng.random() < e:
                b = rng.choice([c for c in "ACGT" if c != b])
            seq.append(b)
            q = max(2, min(41, int(round(-10 * math.log10(e))) + rng.randint(-2, 2)))
            qual.append(chr(33 + q))
        made += 1
        yield "".join(seq), "".join(qual)


# ---------------------------------------------------------------- 정리
def trim_quality(seq, qual, q=20, window=4):
    """3' 끝에서부터 창 평균 품질이 q보다 낮은 동안 잘라 냅니다 (fastp의 끝 자르기와 같은 생각)."""
    end = len(seq)
    while end >= window:
        avg = sum(ord(c) - 33 for c in qual[end - window:end]) / window
        if avg >= q:
            break
        end -= 1
    return seq[:end], qual[:end]


def trim_adapter(seq, qual, adapter=ADAPTER, min_overlap=5):
    """리드 안에서 어댑터가 시작하는 곳을 찾아 그 앞에서 자릅니다. 리드 끝에 걸친 짧은 조각도 봅니다."""
    i = seq.find(adapter)
    if i >= 0:
        return seq[:i], qual[:i]
    for n in range(len(adapter) - 1, min_overlap - 1, -1):    # 끝에 어댑터 앞부분만 걸친 경우
        if seq.endswith(adapter[:n]):
            return seq[:-n], qual[:-n]
    return seq, qual


# ---------------------------------------------------------------- 버블
def walk_branch(start, kset, max_steps):
    """갈림길 하나에서 한 갈래를 따라 합류 지점(들어오는 화살표가 2개 이상인 k-mer) 직전까지 걷습니다."""
    path = [start]
    cur = start
    for _ in range(max_steps):
        nxt = successors(cur, kset)
        if len(nxt) != 1:
            return path, None
        n = nxt[0]
        if len(predecessors(n, kset)) >= 2:
            return path, n                                   # 합류 지점
        path.append(n)
        cur = n
    return path, None


def pop_bubbles(kset, counts, k):
    """두 갈래가 길이 2k 안에서 같은 k-mer로 다시 만나면 버블. 평균 횟수가 낮은 갈래의 k-mer를 지웁니다.
    두 벌 유전체의 이형접합 SNP가 만든 버블은 이렇게 한쪽 벌로 합쳐집니다(콘티그가 하나로 이어짐)."""
    popped = 0
    for x in list(kset):
        if x not in kset:
            continue
        for fw in (x, revcomp(x)):
            nxt = successors(fw, kset)
            if len(nxt) != 2:
                continue
            a, ea = walk_branch(nxt[0], kset, 2 * k)
            b, eb = walk_branch(nxt[1], kset, 2 * k)
            if ea is None or ea != eb:
                continue
            ca = sum(counts.get(canonical(y), 0) for y in a) / len(a)
            cb = sum(counts.get(canonical(y), 0) for y in b) / len(b)
            weak = a if ca < cb else b
            for y in weak:
                kset.discard(canonical(y))
            popped += 1
            break
    return popped
