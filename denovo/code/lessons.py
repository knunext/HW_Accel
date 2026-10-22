"""lessons.py — 2~5장에서 손으로 따라 해 보는 작은 함수들. 페이지의 코드 상자와 같은 코드입니다."""
import collections
import math


# ---------------------------------------------------------------- 1장. 커버리지
def coverage_stats(G, N, L):
    """Lander–Waterman 어림: 평균 커버리지 c, 안 읽힌 염기 비율, 예상 섬(콘티그) 개수."""
    c = N * L / G
    uncovered = math.exp(-c)          # 한 염기가 한 번도 안 덮일 확률
    islands = N * math.exp(-c)        # 리드 하나가 섬의 오른쪽 끝일 확률 × N
    return c, uncovered, islands


# ---------------------------------------------------------------- 2장. OLC
def overlap(a, b, min_len=3):
    """a의 끝(접미사)과 b의 앞(접두사)이 겹치는 가장 긴 길이. 없으면 0."""
    for n in range(min(len(a), len(b)) - 1, min_len - 1, -1):
        if a.endswith(b[:n]):
            return n
    return 0


def greedy_assemble(reads, min_len=3):
    """가장 길게 겹치는 쌍부터 합치기를 더 합칠 쌍이 없을 때까지 되풀이합니다."""
    reads = list(dict.fromkeys(reads))                 # 중복 제거, 순서 유지
    reads = [r for r in reads if not any(r != s and r in s for s in reads)]
    while True:
        best = (0, None, None)
        for i, a in enumerate(reads):                  # 모든 쌍을 비교: n(n-1)번
            for j, b in enumerate(reads):
                if i != j:
                    n = overlap(a, b, min_len)
                    if n > best[0]:
                        best = (n, i, j)
        n, i, j = best
        if n == 0:
            return reads                               # 남은 것들이 콘티그
        merged = reads[i] + reads[j][n:]
        reads = [r for t, r in enumerate(reads) if t not in (i, j)] + [merged]


# ---------------------------------------------------------------- 3장. de Bruijn 그래프
def kmers(text, k):
    """길이 k 창을 한 칸씩 밀며 잘라 낸 조각들."""
    return [text[i:i + k] for i in range(len(text) - k + 1)]


def de_bruijn(text, k):
    """k-mer 하나 = 간선 하나: 앞 k-1 글자(노드) → 뒤 k-1 글자(노드).
    같은 (k-1)-mer는 같은 노드가 되므로 저절로 '접착'됩니다."""
    graph = collections.defaultdict(list)
    for km in kmers(text, k):
        graph[km[:-1]].append(km[1:])
    return graph


# ---------------------------------------------------------------- 4장. 오일러 경로
def find_start(graph):
    """나가는 간선이 들어오는 간선보다 하나 많은 노드에서 출발합니다."""
    out_deg = {u: len(vs) for u, vs in graph.items()}
    in_deg = collections.Counter(v for vs in graph.values() for v in vs)
    for u in graph:
        if out_deg[u] - in_deg[u] == 1:
            return u
    return next(iter(graph))                           # 모두 균형이면 아무 데서나(순환)


def euler_path(graph):
    """Hierholzer 알고리즘: 막힐 때까지 걷고, 막히면 되돌아오며 경로에 적습니다."""
    g = {u: list(vs) for u, vs in graph.items()}       # 쓴 간선을 지우기 위해 복사
    stack = [find_start(g)]
    path = []
    while stack:
        u = stack[-1]
        if g.get(u):                                   # 아직 안 쓴 간선이 있으면
            stack.append(g[u].pop(0))                  #   그 간선을 타고 전진
        else:                                          # 막혔으면
            path.append(stack.pop())                   #   경로에 적고 한 칸 후퇴
    return path[::-1]


def spell(path):
    """노드 경로를 서열로: 첫 노드 + 이후 노드의 마지막 글자."""
    return path[0] + "".join(v[-1] for v in path[1:])


if __name__ == "__main__":
    text = "TAATGCCATGGGATGTT"
    g = de_bruijn(text, 3)
    p = euler_path(g)
    print(dict(g))
    print(p, spell(p), spell(p) == text)
    reads = ["ATGCGT", "GCGTAC", "GTACCA", "ACCATG"]
    print(greedy_assemble(reads))
    print(coverage_stats(5386, 1077, 150))


# ---------------------------------------------------------------- 8장. DPU 나누기 흉내
def mix64(x):
    """64비트 섞기 (splitmix64 마무리). 비슷한 k-mer도 전혀 다른 값이 됩니다."""
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9 & 0xFFFFFFFFFFFFFFFF
    x = (x ^ (x >> 27)) * 0x94D049BB133111EB & 0xFFFFFFFFFFFFFFFF
    return x ^ (x >> 31)


def enc2(s):
    """A=0 C=1 G=2 T=3으로 읽은 4진수 = 2비트씩 붙인 정수."""
    return int(s.translate(str.maketrans("ACGT", "0123")), 4)


def revcomp_s(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def owner_dpu(kmer, n_dpu=256, how="hash", m=11):
    """정규형 k-mer를 맡을 DPU 번호. how: 'prefix' | 'hash' | 'minimizer'"""
    if how == "prefix":                              # 앞 4글자 = 256구역
        return enc2(kmer[:4])
    if how == "hash":                                # k-mer 전체의 해시
        return mix64(enc2(kmer)) % n_dpu
    mmers = (kmer[i:i + m] for i in range(len(kmer) - m + 1))
    best = min(mix64(enc2(min(s, revcomp_s(s)))) for s in mmers)
    return best % n_dpu                              # 해시가 가장 작은 m-mer(미니마이저)


def dpu_loads(kmers, n_dpu=256, how="hash"):
    """DPU마다 맡은 k-mer 수와, 가장 바쁜 DPU가 평균의 몇 배인지."""
    load = [0] * n_dpu
    for km in kmers:
        load[owner_dpu(km, n_dpu, how)] += 1
    return load, max(load) / (sum(load) / n_dpu)


# ---------------------------------------------------------------- 9장. 슈퍼-k-mer
def superkmers(read, k=31, m=11):
    """미니마이저가 같은 연속 k-mer를 한 덩어리로 묶습니다.
    덩어리 하나 = (미니마이저 해시, 서열). 서열 길이 = k - 1 + 묶인 k-mer 수."""
    L = len(read)
    rc = revcomp_s(read)
    mh = []                                            # 위치 j의 정규형 m-mer 해시
    for j in range(L - m + 1):
        fw, bw = read[j:j + m], rc[L - m - j:L - j]
        mh.append(mix64(enc2(fw if fw <= bw else bw)))
    out, start, cur = [], 0, None
    for i in range(L - k + 1):
        mn = min(mh[i:i + k - m + 1])                  # k-mer i의 미니마이저
        if cur is not None and mn != cur:
            out.append((cur, read[start:i - 1 + k]))   # k-mer start..i-1을 담은 서열
            start = i
        cur = mn
    out.append((cur, read[start:L]))
    return out
