// count_pe.cpp — Alveo U50 의사 채널 하나를 맡는 k-mer 세기 엔진 (Vitis HLS 스케치)
// 32개를 복제해 의사 채널 0..31에 하나씩 잇고, 앞단 분배기가 hash(x) % 32로 k-mer를 나눠 준다.
#include <ap_int.h>
#include <hls_stream.h>

#define LOG_SLOTS 23          // 256 MB / 32 B = 2^23 칸 (칸 하나 = AXI 256비트 한 번)
#define FWD 8                 // 최근에 쓴 칸을 들고 있는 수: 읽기→쓰기 지연(반복 횟수)보다 커야 한다

typedef ap_uint<256> bucket_t;            // [63:0] key, [95:64] count, 나머지 비움

static ap_uint<LOG_SLOTS> slot_of(ap_uint<64> x) {
    ap_uint<64> h = x ^ (x >> 29);
    h = h * ap_uint<64>("0xBF58476D1CE4E5B9", 16);        // FPGA는 DSP 곱셈기가 넉넉하다
    return h(63, 64 - LOG_SLOTS);
}

void count_pe(hls::stream<ap_uint<64> > &in, bucket_t *hbm, unsigned n) {
#pragma HLS INTERFACE m_axi port=hbm bundle=pc offset=slave num_read_outstanding=32 num_write_outstanding=32
    ap_uint<LOG_SLOTS> fa[FWD];           // 포워딩 레지스터: 최근에 쓴 칸 주소와 값
    bucket_t fv[FWD];
    bool ok[FWD];
#pragma HLS ARRAY_PARTITION variable=fa complete
#pragma HLS ARRAY_PARTITION variable=fv complete
#pragma HLS ARRAY_PARTITION variable=ok complete
    for (int j = 0; j < FWD; j++) ok[j] = false;

    for (unsigned i = 0; i < n; i++) {
#pragma HLS PIPELINE II=1
#pragma HLS DEPENDENCE variable=hbm inter false   // 메모리를 거친 의존은 아래 포워딩이 책임진다
        ap_uint<64> x = in.read();
        ap_uint<LOG_SLOTS> a = slot_of(x);
        bucket_t b = hbm[a];                          // HBM에서 옛 값 (아직 갱신 전일 수 있다)
        for (int j = FWD - 1; j >= 0; j--)            // 더 최근 것(j가 작은 것)이 마지막에 이긴다
            if (ok[j] && fa[j] == a) b = fv[j];
        ap_uint<64> key = b(63, 0);
        ap_uint<32> cnt = b(95, 64);
        if (cnt == 0 || key == x) { b(63, 0) = x; b(95, 64) = cnt + 1; }
        // 생략: 다른 k-mer가 이미 앉은 칸(충돌) — 칸 하나에 key 두 개를 담거나, 넘침 큐로 보내 다음 회차에 처리
        hbm[a] = b;
        for (int j = FWD - 1; j > 0; j--) { fa[j] = fa[j - 1]; fv[j] = fv[j - 1]; ok[j] = ok[j - 1]; }
        fa[0] = a; fv[0] = b; ok[0] = true;
    }
}
