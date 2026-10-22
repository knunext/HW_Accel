"""bucket_compact.py — 구역(버킷)마다 따로 압축하고 마지막에 이어 붙여 유니티그를 만듭니다 (14장)

GPU·FPGA 구현의 기준 코드입니다. 구역 하나가 장치의 작업 하나예요.

  1. 나누기   : k-mer를 앞 (k-1)-mer와 뒤 (k-1)-mer의 미니마이저 구역에 둘 다 넣는다 (최대 두 곳)
  2. 구역 안  : 그 구역이 '주인'인 (k-1)-mer마다 왼쪽·오른쪽으로 이어지는 k-mer가 하나씩뿐인지 본다
                → 그렇다면 두 k-mer는 같은 유니티그 (갈림길 없는 이음)
  3. 잇기     : 모든 구역이 내놓은 이음을 union-find로 합친다 (CPU)

(k-1)-mer s로 이어지는 k-mer는 모두 s를 앞이나 뒤에 품고 있고, 그래서 모두 s의 구역에 들어 있습니다.
그러니 2단계는 다른 구역을 볼 필요가 없어요. 결과는 mini_assembler.build_unitigs와 같은 k-mer 묶음입니다.
"""
import sys

from mini_assembler import canonical, mix64

TR = str.maketrans("ACGT", "0123")


def minimizer_hash(s, m=11):
    """s 안의 정규형 m-mer 해시 중 최솟값. 가닥과 무관합니다."""
    best = None
    for i in range(len(s) - m + 1):
        h = mix64(int(canonical(s[i:i + m]).translate(TR), 4))
        if best is None or h < best:
            best = h
    return best


def side_bucket(s, n_bins, m=11):
    return minimizer_hash(s, m) % n_bins


def partition(kset, n_bins, m=11):
    """1. 나누기: 앞뒤 (k-1)-mer의 구역에 넣습니다. 두 구역이 같으면 한 번만."""
    buckets = [set() for _ in range(n_bins)]
    for x in kset:
        buckets[side_bucket(x[:-1], n_bins, m)].add(x)
        buckets[side_bucket(x[1:], n_bins, m)].add(x)
    return buckets


def bucket_links(bucket, b, n_bins, m=11):
    """2. 구역 안: 이 구역이 주인인 (k-1)-mer 중 갈림길 없는 것마다 (왼쪽 k-mer, 오른쪽 k-mer)."""
    seen = set()
    for x in bucket:
        for s in (x[:-1], x[1:]):
            cs = canonical(s)
            if cs in seen:
                continue
            seen.add(cs)
            if side_bucket(cs, n_bins, m) != b:      # 주인 구역이 아닌 쪽은 건너뜀
                continue
            left = [y for y in (canonical(a + cs) for a in "ACGT") if y in bucket]
            right = [y for y in (canonical(cs + c) for c in "ACGT") if y in bucket]
            if len(left) == 1 and len(right) == 1:
                yield left[0], right[0]


class UnionFind:
    def __init__(self):
        self.p = {}

    def find(self, x):
        p = self.p
        p.setdefault(x, x)
        while p[x] != x:
            p[x] = p[p[x]]                           # 경로 반으로 줄이기
            x = p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def bucket_unitigs(kset, n_bins=64, m=11):
    """3. 잇기까지. 돌려주는 값: 유니티그마다 k-mer 묶음, 그리고 통계."""
    buckets = partition(kset, n_bins, m)
    uf = UnionFind()
    links = 0
    for b, bk in enumerate(buckets):                 # 구역마다 독립: GPU 블록이나 장치 하나에 맡길 수 있음
        for x, y in bucket_links(bk, b, n_bins, m):
            uf.union(x, y)
            links += 1
    groups = {}
    for x in kset:
        groups.setdefault(uf.find(x), set()).add(x)
    stats = dict(kmers=len(kset), placed=sum(len(bk) for bk in buckets),
                 largest=max(len(bk) for bk in buckets), links=links, unitigs=len(groups))
    return list(groups.values()), stats


if __name__ == "__main__":
    import mini_assembler as ma
    reads = sys.argv[1]
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 31
    counts = ma.count_kmers(ma.read_fastq(reads), k)
    thr = ma.auto_threshold(ma.kmer_spectrum(counts))
    kset = {x for x, c in counts.items() if c >= thr}
    groups, st = bucket_unitigs(kset, 64, m=min(11, k - 4))
    ref = {frozenset(canonical(u[i:i + k]) for i in range(len(u) - k + 1)) for u in ma.build_unitigs(set(kset), k)}
    print(st, "중복 배치", round(st["placed"] / st["kmers"], 3), "전역 판과 같음:", ref == {frozenset(g) for g in groups})
