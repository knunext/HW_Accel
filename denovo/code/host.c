// host.c — k-mer를 DPU별로 나눠 보내고, 세게 하고, 결과를 받는다 (UPMEM SDK 스케치)
#include <dpu.h>
#include <stdint.h>
#include <stdlib.h>

#define NR_DPUS 256                    // 우리 서버: PIM 16 GB = DPU 256개
#define DPU_BINARY "./kmer_dpu"

static inline uint64_t mix64(uint64_t x) {   // 호스트는 곱셈이 싸니 splitmix64를 쓴다
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

int main(void) {
    struct dpu_set_t set, dpu;
    uint32_t i;
    DPU_ASSERT(dpu_alloc(NR_DPUS, NULL, &set));
    DPU_ASSERT(dpu_load(set, DPU_BINARY, NULL));

    // 1) 리드에서 정규형 k-mer를 2비트 정수로 꺼내 owner = mix64(x) % NR_DPUS 묶음에 담는다
    uint64_t *buf[NR_DPUS];
    uint32_t n_kmers[NR_DPUS], max_n = 0;
    /* ... fill_buckets(reads, buf, n_kmers) ... */
    for (i = 0; i < NR_DPUS; i++) if (n_kmers[i] > max_n) max_n = n_kmers[i];

    // 2) 개수와 k-mer 묶음을 병렬로 보낸다. 병렬 전송은 모든 DPU가 같은 크기 → max_n에 맞춘다
    DPU_FOREACH(set, dpu, i) { DPU_ASSERT(dpu_prepare_xfer(dpu, &n_kmers[i])); }
    DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "n_kmers", 0, sizeof(uint32_t), DPU_XFER_DEFAULT));
    DPU_FOREACH(set, dpu, i) { DPU_ASSERT(dpu_prepare_xfer(dpu, buf[i])); }
    DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "kmers", 0, max_n * sizeof(uint64_t), DPU_XFER_DEFAULT));

    // 3) 모든 DPU가 동시에 센다
    DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));

    // 4) 문턱을 넘긴 k-mer 개수부터 받고, 그만큼만 가져온다
    uint32_t n_solid[NR_DPUS];
    DPU_FOREACH(set, dpu, i) { DPU_ASSERT(dpu_prepare_xfer(dpu, &n_solid[i])); }
    DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "n_solid", 0, sizeof(uint32_t), DPU_XFER_DEFAULT));
    /* ... DPU마다 dpu_copy_from(dpu, "solid", 0, out[i], n_solid[i] * 16) ... */

    DPU_ASSERT(dpu_free(set));
    return 0;
}
