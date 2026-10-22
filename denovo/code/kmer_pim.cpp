// kmer_pim.cpp — 적층 메모리 x86 PIM 흉내: 스레드 하나 = 볼트 하나 (Sniper 실험용 스케치)
// 빌드: g++ -O2 -pthread -I$SNIPER_ROOT/include kmer_pim.cpp -o kmer_pim
// 실행 예: run-sniper -n 32 -c <PIM 설정> --roi -- ./kmer_pim kmers.bin message
#include <atomic>
#include <cstdint>
#include <cstring>
#include <thread>
#include <vector>
#include "sim_api.h"                         // SimRoiStart(), SimRoiEnd()

constexpr int V = 32;                        // 볼트 수 = 스레드 수
constexpr uint32_t SLOTS = 1u << 20;         // 볼트마다 표 칸 수 (2의 거듭제곱)

struct Slot { std::atomic<uint64_t> key{0}; std::atomic<uint32_t> count{0}; };
struct alignas(64) Vault {
    Slot table[SLOTS];                       // 이 볼트의 메모리에 놓인다고 가정 (첫 접촉으로 배치)
    std::vector<uint64_t> out[V];            // 메시지 방식: 볼트 d로 보낼 k-mer 묶음
};
static Vault *vault[V];

static inline uint64_t mix64(uint64_t x) {
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
static inline int owner(uint64_t x) { return mix64(x) % V; }   // 걷기까지 보려면 미니마이저 해시로 바꾼다

static void insert(Vault *v, uint64_t x) {                  // 선형 탐사 + CAS (key 0은 빈 칸 표시)
    uint32_t p = (mix64(x) >> 8) & (SLOTS - 1);
    for (;;) {
        uint64_t k = v->table[p].key.load(std::memory_order_relaxed);
        if (k == 0 && v->table[p].key.compare_exchange_strong(k, x)) k = x;   // 실패하면 k에 남이 넣은 값
        if (k == x) { v->table[p].count.fetch_add(1, std::memory_order_relaxed); return; }
        p = (p + 1) & (SLOTS - 1);
    }
}

// 공유 메모리 방식: 아무 스레드나 주인 볼트의 표를 직접 건드린다 (대부분 원격 접근)
static void count_shared(int t, const uint64_t *km, size_t n) {
    for (size_t i = t; i < n; i += V) insert(vault[owner(km[i])], km[i]);
}

// 메시지 방식: (1) 내 몫을 주인별 묶음에 담고 (2) 장벽 뒤 내게 온 묶음만 내 표에 센다 (모두 로컬 접근)
static std::atomic<int> arrived{0};
static void count_message(int t, const uint64_t *km, size_t n) {
    for (size_t i = t; i < n; i += V) vault[t]->out[owner(km[i])].push_back(km[i]);
    arrived.fetch_add(1);
    while (arrived.load() < V) { }                           // 간단한 장벽
    for (int s = 0; s < V; s++)
        for (uint64_t x : vault[s]->out[t]) insert(vault[t], x);
}

int main(int argc, char **argv) {
    std::vector<uint64_t> km;                                // ... argv[1]에서 2비트 정규형 k-mer를 읽어 채운다 ...
    bool message = argc > 2 && !strcmp(argv[2], "message");
    for (int v = 0; v < V; v++) vault[v] = new Vault();     // 실제로는 볼트 v의 코어가 직접 할당해 첫 접촉 배치
    SimRoiStart();
    std::vector<std::thread> th;
    for (int t = 0; t < V; t++)
        th.emplace_back(message ? count_message : count_shared, t, km.data(), km.size());
    for (auto &x : th) x.join();
    SimRoiEnd();
    return 0;
}
