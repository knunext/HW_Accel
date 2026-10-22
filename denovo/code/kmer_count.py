"""kmer_count.py — k-mer 세기를 빠르게 만드는 방법들 (8장)

  1. 2비트로 굴리기      : 문자열 대신 64비트 정수, 한 글자에 시프트 한 번
  2. 사전(해시 표)       : k-mer마다 표의 아무 자리를 한 번 건드림
  3. 정렬 후 세기        : 모두 모아 정렬하고, 같은 값이 이어진 길이를 셈
  4. 기수 정렬(LSD)      : 비교 없이 16비트씩 네 번 나눠 담기
  5. 나누고 세기         : 미니마이저로 구역을 나눠, 구역마다 따로 정렬
  6. 스펙트럼으로 유전체 크기 어림

2~5는 같은 답을 냅니다. numpy가 필요한 함수는 3~5입니다.
"""
import collections

ENC = {"A": 0, "C": 1, "G": 2, "T": 3}


# ---------------------------------------------------------------- 1. 2비트로 굴리기
def kmers_2bit(read, k):
    """정규형 k-mer를 정수로 하나씩 내놓습니다 (k <= 31).
    fw: 왼쪽으로 2비트 밀고 새 염기를 오른쪽에 넣음
    rc: 오른쪽으로 2비트 밀고 새 염기의 상보(3 - b)를 왼쪽 끝에 넣음"""
    mask = (1 << (2 * k)) - 1
    shift = 2 * (k - 1)
    fw = rc = 0
    n = 0                                   # 지금까지 이어 읽은 염기 수 (N을 만나면 0)
    for ch in read:
        b = ENC.get(ch)
        if b is None:
            n = 0
            continue
        fw = ((fw << 2) | b) & mask
        rc = (rc >> 2) | ((3 - b) << shift)
        n += 1
        if n >= k:
            yield fw if fw < rc else rc     # 정규형 = 둘 중 작은 정수


def count_dict(reads, k):
    """해시 표(파이썬 사전)로 세기."""
    counts = collections.Counter()
    for r in reads:
        for x in kmers_2bit(r, k):
            counts[x] += 1
    return counts


# ---------------------------------------------------------------- 3~5. numpy로 한꺼번에
def encode_reads(reads, k):
    """같은 길이 리드들을 2차원 배열로 바꿔 모든 정규형 k-mer를 한 번에 계산합니다."""
    import numpy as np
    lut = np.full(256, 255, dtype=np.uint8)
    for ch, v in ENC.items():
        lut[ord(ch)] = v
    B = lut[np.frombuffer("".join(reads).encode(), dtype=np.uint8)].reshape(len(reads), -1).astype(np.uint64)
    L = B.shape[1]
    W = L - k + 1
    fw = np.zeros((len(reads), W), dtype=np.uint64)
    rc = np.zeros((len(reads), W), dtype=np.uint64)
    for j in range(k):                                   # k번의 시프트로 모든 위치를 동시에
        fw = (fw << np.uint64(2)) | B[:, j:j + W]
        rc = rc | ((np.uint64(3) - B[:, j:j + W]) << np.uint64(2 * j))
    return np.minimum(fw, rc).ravel()


def count_runs(sorted_keys):
    """정렬된 배열에서 같은 값이 이어진 길이 = 횟수."""
    import numpy as np
    if sorted_keys.size == 0:
        return sorted_keys, sorted_keys
    edge = np.flatnonzero(sorted_keys[1:] != sorted_keys[:-1]) + 1
    starts = np.concatenate(([0], edge))
    counts = np.diff(np.concatenate((starts, [sorted_keys.size])))
    return sorted_keys[starts], counts


def count_sort(keys):
    """정렬 후 세기 (numpy 기본 정렬)."""
    import numpy as np
    return count_runs(np.sort(keys))


def radix_sort_lsd(keys, bits=16, key_bits=62):
    """LSD 기수 정렬: 낮은 자리 16비트부터 안정적으로 나눠 담기를 되풀이합니다.
    비교를 하지 않고, 매번 배열을 앞에서 뒤로 순서대로 읽고 씁니다."""
    import numpy as np
    out = keys
    for shift in range(0, key_bits, bits):
        digit = ((out >> np.uint64(shift)) & np.uint64((1 << bits) - 1)).astype(np.uint16)
        order = np.argsort(digit, kind="stable")        # 16비트 정수의 안정 정렬 = 계수(counting) 정렬
        out = out[order]
    return out


def count_radix(keys, k=31):
    return count_runs(radix_sort_lsd(keys, 16, 2 * k))


def minimizer_bins(reads, k, m=11, n_bins=64):
    """k-mer마다 미니마이저(해시가 가장 작은 정규형 m-mer)로 구역 번호를 정합니다."""
    import numpy as np
    mm = encode_reads(reads, m).reshape(len(reads), -1)               # 위치별 정규형 m-mer
    h = (mm * np.uint64(0x9E3779B97F4A7C15)) >> np.uint64(40)       # 곱셈 해시 (넘침은 버림)
    W = len(reads[0]) - k + 1
    best = h[:, 0:W].copy()
    for j in range(1, k - m + 1):                                     # 창 안의 최솟값
        best = np.minimum(best, h[:, j:j + W])
    return (best % np.uint64(n_bins)).astype(np.int64).ravel()


def count_partitioned(reads, k=31, n_bins=64):
    """나누고 세기: 구역마다 따로 기수 정렬. 한 번에 다루는 배열이 1/n_bins로 작아집니다."""
    import numpy as np
    keys = encode_reads(reads, k)
    bins = minimizer_bins(reads, k, n_bins=n_bins)
    order = np.argsort(bins, kind="stable")
    keys, bins = keys[order], bins[order]
    cuts = np.searchsorted(bins, np.arange(n_bins + 1))
    outs = [count_radix(keys[cuts[i]:cuts[i + 1]], k) for i in range(n_bins)]
    largest = int(np.max(np.diff(cuts)))
    return outs, largest


# ---------------------------------------------------------------- 6. 스펙트럼으로 유전체 크기 어림
def genome_size(spectrum, threshold):
    """오류 봉우리를 뺀 k-mer 총수 / 진짜 봉우리 높이(가장 흔한 횟수) ≈ 유전체 길이."""
    peak = max(range(threshold, len(spectrum)), key=lambda c: spectrum[c])
    total = sum(c * spectrum[c] for c in range(threshold, len(spectrum)))
    return total / peak, peak
