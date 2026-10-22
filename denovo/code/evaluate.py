"""evaluate.py — N50 너머: 콘티그를 참조에 맞춰 보고 NGA50과 오조립을 셉니다 (16장)

QUAST가 하는 일을 아주 작게 흉내 낸 것입니다.
  1. 참조에서 한 번만 나오는 k-mer를 닻(anchor)으로 콘티그의 위치를 찾고
  2. 콘티그 위치 - 참조 위치(대각선)가 같은 닻끼리 묶어 '맞는 덩어리'를 만들고
  3. 대각선이나 가닥이 바뀌는 곳을 오조립으로 셉니다.
"""
import collections

from mini_assembler import canonical, n50, revcomp, successors


def unique_kmer_index(ref, k):
    """참조에서 딱 한 번 나오는 정규형 k-mer → (위치, 가닥)."""
    seen = collections.Counter()
    pos = {}
    for i in range(len(ref) - k + 1):
        km = ref[i:i + k]
        c = canonical(km)
        seen[c] += 1
        pos[c] = (i, km == c)
    return {c: p for c, p in pos.items() if seen[c] == 1}


def align_blocks(contig, index, k, max_gap=1000):
    """콘티그를 '맞는 덩어리' 목록으로: (콘티그 시작, 끝, 참조 시작, 가닥)."""
    anchors = []
    for i in range(len(contig) - k + 1):
        km = contig[i:i + k]
        c = canonical(km)
        if c in index:
            rpos, ref_is_canon = index[c]
            strand = (km == c) == ref_is_canon               # 참조와 같은 가닥이면 True
            diag = rpos - i if strand else rpos + i
            anchors.append((i, rpos, strand, diag))
    blocks = []
    for a in anchors:
        if blocks and blocks[-1][3] == a[2] and abs(blocks[-1][4] - a[3]) <= 2 and a[0] - blocks[-1][1] <= max_gap:
            b = blocks[-1]
            blocks[-1] = (b[0], a[0], min(b[2], a[1]), b[3], b[4])
        else:
            blocks.append((a[0], a[0], a[1], a[2], a[3]))
    return [(b[0], b[1] + k, b[2], b[3]) for b in blocks]


def evaluate_contigs(contigs, ref, k=31):
    """N50, NG50, NGA50, 오조립 수, 참조를 덮은 비율."""
    idx = unique_kmer_index(ref, k)
    lens = [len(c) for c in contigs]
    blocks, mis, covered = [], 0, bytearray(len(ref))
    for c in contigs:
        bl = align_blocks(c, idx, k)
        mis += max(0, len(bl) - 1)                          # 덩어리가 둘 이상 = 이음매가 틀린 곳
        for s, e, r, _ in bl:
            blocks.append(e - s)
            covered[r:r + (e - s)] = b"\x01" * min(e - s, len(ref) - r)
    G = len(ref)
    def ngx(L):
        acc = 0
        for x in sorted(L, reverse=True):
            acc += x
            if acc * 2 >= G:
                return x
        return 0
    return dict(contigs=len(contigs), n50=n50(lens), ng50=ngx(lens), nga50=ngx(blocks),
                misassemblies=mis, genome_fraction=sum(covered) / G)


def greedy_extend(start, kset, counts, used):
    """갈림길에서도 멈추지 않고, 아직 안 쓴 쪽 중 가장 자주 나온 k-mer로 계속 걷습니다."""
    seq, cur = start, start
    while True:
        nxt = [n for n in successors(cur, kset) if canonical(n) not in used]
        if not nxt:
            return seq
        n = max(nxt, key=lambda y: (counts.get(canonical(y), 0), y))
        used.add(canonical(n))
        seq += n[-1]
        cur = n


def greedy_contigs(kset, counts, k):
    """비교용 '과감한' 조립기: 씨앗에서 양쪽으로 갈림길을 무시하고 걷습니다.
    콘티그는 길어지지만 반복 앞뒤를 틀리게 이을 수 있어요."""
    used, out = set(), []
    for seed in sorted(kset):
        if seed in used:
            continue
        used.add(seed)
        right = greedy_extend(seed, kset, counts, used)
        left = greedy_extend(revcomp(seed), kset, counts, used)
        out.append(revcomp(left)[:-k] + right)
    return out
