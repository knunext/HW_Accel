// kmer_dpu.c — DPU 하나가 받은 k-mer를 MRAM 해시 표에 센다 (스케치)
// 컴파일: dpu-upmem-dpurte-clang -DNR_TASKLETS=16 -O2 -o kmer_dpu kmer_dpu.c
#include <stdint.h>
#include <defs.h>      // me()
#include <mram.h>      // mram_read, mram_write (8바이트 배수, 최대 2048바이트)
#include <barrier.h>

#define MAX_KMERS   (3u << 20)     // 받을 k-mer 최대 300만 개 = 24 MB
#define SLOTS       (1u << 20)     // 표 100만 칸 × 16 B = 16 MB
#define SPAN        (SLOTS / NR_TASKLETS)
#define BLOCK       256            // 한 번에 WRAM으로 읽는 k-mer 수 = 2 KB

typedef struct { uint64_t key; uint64_t count; } slot_t;

__host uint32_t n_kmers;
__host uint32_t n_solid;
__mram_noinit uint64_t kmers[MAX_KMERS];
__mram_noinit slot_t table[SLOTS];
__dma_aligned uint64_t cache[NR_TASKLETS][BLOCK];
__dma_aligned slot_t cur[NR_TASKLETS];
BARRIER_INIT(sync, NR_TASKLETS);

static inline uint32_t xhash(uint64_t x) {   // 곱셈 없이 시프트와 XOR만 (DPU는 곱셈이 느리다)
    x ^= x >> 31; x ^= x << 17; x ^= x >> 23;
    return (uint32_t)x ^ (uint32_t)(x >> 32);
}

int main(void) {
    uint32_t t = me(), lo = t * SPAN;

    // 0) 내 구역을 0으로 비운다 (__mram_noinit는 쓰레기 값으로 시작)
    cur[t].key = 0; cur[t].count = 0;
    for (uint32_t p = lo; p < lo + SPAN; p++) mram_write(&cur[t], &table[p], sizeof(slot_t));
    barrier_wait(&sync);

    // 1) 모든 tasklet이 같은 블록을 읽고, 해시가 내 구역에 떨어지는 k-mer만 센다 → 잠금 불필요
    for (uint32_t base = 0; base < n_kmers; base += BLOCK) {
        uint32_t n = (n_kmers - base < BLOCK) ? n_kmers - base : BLOCK;
        mram_read(&kmers[base], cache[t], BLOCK * sizeof(uint64_t));
        for (uint32_t j = 0; j < n; j++) {
            uint64_t x = cache[t][j];
            uint32_t hv = xhash(x);
            if ((hv & (NR_TASKLETS - 1)) != t) continue;
            uint32_t p = lo + ((hv >> 4) & (SPAN - 1));            // NR_TASKLETS = 16 = 2^4
            for (;;) {                                              // 선형 탐사
                mram_read(&table[p], &cur[t], sizeof(slot_t));
                if (cur[t].count == 0 || cur[t].key == x) {
                    cur[t].key = x; cur[t].count++;
                    mram_write(&cur[t], &table[p], sizeof(slot_t));
                    break;
                }
                p = (p + 1 < lo + SPAN) ? p + 1 : lo;
            }
        }
    }
    barrier_wait(&sync);

    // 2) 문턱(예: 4) 이상인 칸만 앞으로 모아 n_solid에 적는다 (생략: tasklet별 개수 → 접두합 → 복사)
    return 0;
}
