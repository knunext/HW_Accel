// count_bench.cpp — 같은 정규형 31-mer 2,400만 개를 C++로 네 가지 방법으로 센다 (1코어)
// g++ -O3 -march=native count_bench.cpp -o count_bench && ./count_bench keys.bin bins.bin
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>
using namespace std;
using u64 = uint64_t;
static double now() { return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count(); }
struct Res { size_t distinct; u64 chk; };
static Res runs(const u64 *a, size_t n) {                 // 정렬된 배열의 덩어리 세기
    Res r{0, 0};
    for (size_t i = 0; i < n;) { size_t j = i; while (j < n && a[j] == a[i]) j++; r.distinct++; r.chk += a[i] * (j - i); i = j; }
    return r;
}
static inline u64 mix(u64 x) { x ^= x >> 31; x *= 0x7fb5d329728ea185ULL; x ^= x >> 27; x *= 0x81dadef4bc2dd44dULL; return x ^ (x >> 33); }

// 1) 해시 표: 열린 주소, 선형 탐사, 칸 16바이트
static Res hash_count(const vector<u64> &k, int log_slots) {
    size_t S = size_t(1) << log_slots, M = S - 1;
    vector<u64> key(S, ~0ULL); vector<uint32_t> cnt(S, 0);
    for (u64 x : k) { size_t p = mix(x) & M; while (key[p] != x && key[p] != ~0ULL) p = (p + 1) & M; key[p] = x; cnt[p]++; }
    Res r{0, 0}; for (size_t p = 0; p < S; p++) if (cnt[p]) { r.distinct++; r.chk += key[p] * cnt[p]; }
    return r;
}
// 2) LSD 기수 정렬: bits비트씩 나눠 담기 (계수 정렬을 되풀이)
static void radix(u64 *a, u64 *tmp, size_t n, int bits, int key_bits) {
    size_t B = size_t(1) << bits; vector<size_t> c(B);
    u64 *src = a, *dst = tmp; int passes = 0;
    for (int sh = 0; sh < key_bits; sh += bits, passes++) {
        fill(c.begin(), c.end(), 0);
        for (size_t i = 0; i < n; i++) c[(src[i] >> sh) & (B - 1)]++;
        size_t s = 0; for (size_t b = 0; b < B; b++) { size_t t = c[b]; c[b] = s; s += t; }
        for (size_t i = 0; i < n; i++) dst[c[(src[i] >> sh) & (B - 1)]++] = src[i];
        swap(src, dst);
    }
    if (src != a) memcpy(a, src, n * sizeof(u64));
}
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb"); fseek(f, 0, SEEK_END); size_t n = ftell(f) / 8; rewind(f);
    vector<u64> k(n); fread(k.data(), 8, n, f); fclose(f);
    vector<uint8_t> bin(n); f = fopen(argv[2], "rb"); fread(bin.data(), 1, n, f); fclose(f);
    double t; Res r;
    t = now(); r = hash_count(k, 22); printf("hash %.3f %zu %llu\n", now() - t, r.distinct, (unsigned long long)r.chk);
    { vector<u64> a = k; t = now(); sort(a.begin(), a.end()); r = runs(a.data(), n); printf("stdsort %.3f %zu %llu\n", now() - t, r.distinct, (unsigned long long)r.chk); }
    for (int bits : {8, 11, 16}) { vector<u64> a = k, tmp(n); t = now(); radix(a.data(), tmp.data(), n, bits, 62); r = runs(a.data(), n); printf("radix%d %.3f %zu %llu\n", bits, now() - t, r.distinct, (unsigned long long)r.chk); }
    {   // 4) 64구역으로 나누고 구역마다 기수 정렬(11비트)
        t = now();
        vector<size_t> c(65, 0); for (size_t i = 0; i < n; i++) c[bin[i] + 1]++;
        for (int b = 0; b < 64; b++) c[b + 1] += c[b];
        vector<size_t> pos(c.begin(), c.end() - 1); vector<u64> sc(n);
        for (size_t i = 0; i < n; i++) sc[pos[bin[i]]++] = k[i];
        size_t mx = 0; for (int b = 0; b < 64; b++) mx = max(mx, c[b + 1] - c[b]);
        vector<u64> tmp(mx); Res tot{0, 0};
        for (int b = 0; b < 64; b++) { size_t m = c[b + 1] - c[b]; radix(sc.data() + c[b], tmp.data(), m, 11, 62); Res q = runs(sc.data() + c[b], m); tot.distinct += q.distinct; tot.chk += q.chk; }
        printf("part64 %.3f %zu %llu %zu\n", now() - t, tot.distinct, (unsigned long long)tot.chk, mx);
    }
    {   // 5) 64구역으로 나누고 구역마다 해시 표 (표가 작아 캐시에 가까워짐)
        t = now();
        vector<size_t> c(65, 0); for (size_t i = 0; i < n; i++) c[bin[i] + 1]++;
        for (int b = 0; b < 64; b++) c[b + 1] += c[b];
        vector<size_t> pos(c.begin(), c.end() - 1); vector<u64> sc(n);
        for (size_t i = 0; i < n; i++) sc[pos[bin[i]]++] = k[i];
        Res tot{0, 0};
        for (int b = 0; b < 64; b++) { vector<u64> part(sc.begin() + c[b], sc.begin() + c[b + 1]); Res q = hash_count(part, 16); tot.distinct += q.distinct; tot.chk += q.chk; }
        printf("parthash %.3f %zu %llu\n", now() - t, tot.distinct, (unsigned long long)tot.chk);
    }
    return 0;
}
