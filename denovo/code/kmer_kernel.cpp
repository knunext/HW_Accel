// kmer_kernel.cpp — SmartSSD FPGA: FASTQ 글자 흐름 → 2비트 → k-mer → 미니마이저 → 슈퍼-k-mer (Vitis HLS 스케치)
#include <ap_int.h>
#include <hls_stream.h>
#define K 31
#define M 11
#define W (K - M + 1)              // 한 k-mer 안의 m-mer 수

typedef ap_uint<2 * K> kmer_t;
typedef ap_uint<2 * M> mmer_t;
struct sk_word { ap_uint<64> bits; ap_uint<8> len; ap_uint<1> last; };   // 슈퍼-k-mer 조각

static ap_uint<32> mhash(mmer_t v) {        // FPGA는 곱셈기(DSP)가 넉넉하다
    ap_uint<64> x = v;
    x = (x ^ (x >> 15)) * 0x2C1B3C6DULL;
    return x(31, 0) ^ x(63, 32);
}

void kmer_kernel(hls::stream<ap_uint<8> > &in, hls::stream<sk_word> &out, unsigned n_chars) {
    ap_uint<2> line = 0;                    // FASTQ 줄 번호: 0 이름, 1 서열, 2 '+', 3 품질
    kmer_t fw = 0;                          // k글자 시프트 레지스터
    mmer_t mm = 0;
    ap_uint<32> win[W];                     // 창 안의 m-mer 해시들
#pragma HLS ARRAY_PARTITION variable=win complete
    unsigned nbase = 0;
    ap_uint<32> prev_min = 0xFFFFFFFF;
    sk_word acc = {0, 0, 0};
    for (unsigned i = 0; i < n_chars; i++) {
#pragma HLS PIPELINE II=1
        ap_uint<8> c = in.read();
        if (c == '\n') { line++; if (line == 2) nbase = 0; continue; }   // 서열 줄이 끝나면 레지스터 비움
        if (line != 1) continue;                                         // 서열 줄이 아니면 버린다 (빈 클럭)
        ap_uint<2> b = (c == 'A') ? 0 : (c == 'C') ? 1 : (c == 'G') ? 2 : 3;
        fw = (fw << 2) | b;
        mm = (mm << 2) | b;
        for (int j = 0; j < W - 1; j++) win[j] = win[j + 1];             // 창 밀기 (레지스터 시프트)
        win[W - 1] = mhash(mm);
        nbase++;
        if (nbase < K) continue;
        ap_uint<32> mn = win[0];
        for (int j = 1; j < W; j++) if (win[j] < mn) mn = win[j];        // W개 비교 트리로 합성됨
        if (mn != prev_min && acc.len > 0) { acc.last = 1; out.write(acc); acc.len = 0; }
        prev_min = mn;
        acc.bits = (acc.bits << 2) | b; acc.len++;                       // 실제로는 64비트가 차면 조각을 내보낸다
        // 생략한 것: 새 덩어리는 앞 k-1글자를 다시 담아야 하고, 정규형(역상보 레지스터)과 리드 끝 처리도 필요하다
    }
}
