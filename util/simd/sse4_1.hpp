#ifndef SIMD_SSE4_1_HPP
#define SIMD_SSE4_1_HPP

#include "ssse3.hpp"

#include "smmintrin.h"

namespace simd {

    // load
    static inline m128i movntdqa(m128i * a) {
        return _mm_stream_load_si128(a);
    }

    // pack
    static inline m128i packusdw(m128i a, m128i b) {
        return _mm_packus_epi32(a, b);
    }

    // extract and insert
    template<int N>
    static inline int32_t extractps(m128 a) {
        return _mm_extract_ps(a, N);
    }
    template<int N>
    static inline int pextrb(m128i a) {
        return _mm_extract_epi8(a, N);
    }
    template<int N>
    static inline int32_t pextrd(m128i a) {
        return _mm_extract_epi32(a, N);
    }
    template<int N>
    static inline int64_t pextrq(m128i a) {
        return _mm_extract_epi64(a, N);
    }
    template<int N>
    static inline m128 insertps(m128 a, m128 b) {
        return _mm_insert_ps(a, b, N);
    }
    template<int N>
    static inline m128i pinsrb(m128i a, int b) {
        return _mm_insert_epi8(a, b, N);
    }
    template<int N>
    static inline m128i pinsrd(m128i a, int32_t b) {
        return _mm_insert_epi32(a, b, N);
    }
    template<int N>
    static inline m128i pinsrq(m128i a, int64_t b) {
        return _mm_insert_epi64(a, b, N);
    }

    // blend
    template<int N>
    static inline m128 blendps(m128 a, m128 b) {
        return _mm_blend_ps(a, b, N);
    }
    template<int N>
    static inline m128d blendpd(m128d a, m128d b) {
        return _mm_blend_pd(a, b, N);
    }
    template<int N>
    static inline m128i pblendw(m128i a, m128i b) {
        return _mm_blend_epi16(a, b, N);
    }
    static inline m128 blendvps(m128 a, m128 b, m128 mask) {
        return _mm_blendv_ps(a, b, mask);
    }
    static inline m128d blendvpd(m128d a, m128d b, m128d mask) {
        return _mm_blendv_pd(a, b, mask);
    }
    static inline m128i pblendvb(m128i a, m128i b, m128i mask) {
        return _mm_blendv_epi8(a, b, mask);
    }

    // sign extend
    static inline m128i pmovsxbw(m128i a) {
        return _mm_cvtepi8_epi16(a);
    }
    static inline m128i pmovsxbd(m128i a) {
        return _mm_cvtepi8_epi32(a);
    }
    static inline m128i pmovsxbq(m128i a) {
        return _mm_cvtepi8_epi64(a);
    }
    static inline m128i pmovsxwd(m128i a) {
        return _mm_cvtepi16_epi32(a);
    }
    static inline m128i pmovsxwq(m128i a) {
        return _mm_cvtepi16_epi64(a);
    }
    static inline m128i pmovsxdq(m128i a) {
        return _mm_cvtepi32_epi64(a);
    }

    // zero extend
    static inline m128i pmovzxbw(m128i a) {
        return _mm_cvtepu8_epi16(a);
    }
    static inline m128i pmovzxbd(m128i a) {
        return _mm_cvtepu8_epi32(a);
    }
    static inline m128i pmovzxbq(m128i a) {
        return _mm_cvtepu8_epi64(a);
    }
    static inline m128i pmovzxwd(m128i a) {
        return _mm_cvtepu16_epi32(a);
    }
    static inline m128i pmovzxwq(m128i a) {
        return _mm_cvtepu16_epi64(a);
    }
    static inline m128i pmovzxdq(m128i a) {
        return _mm_cvtepu32_epi64(a);
    }

    // sum of absolute differences
    template<int N>
    static inline m128i mpsadbw(m128i a, m128i b) {
        return _mm_mpsadbw_epu8(a, b, N);
    }

    // multiplication
    static inline m128i pmuldq(m128i a, m128i b) {
        return _mm_mul_epi32(a, b);
    }
    static inline m128i pmulld(m128i a, m128i b) {
        return _mm_mullo_epi32(a, b);
    }

    // dot product
    template<int N>
    static inline m128 dpps(m128 a, m128 b) {
        return _mm_dp_ps(a, b, N);
    }
    template<int N>
    static inline m128d dppd(m128d a, m128d b) {
        return _mm_dp_pd(a, b, N);
    }

    // round
    static inline m128 ceilss(m128 a, m128 b) {
        return _mm_ceil_ss(a, b);
    }
    static inline m128 ceilps(m128 a) {
        return _mm_ceil_ps(a);
    }
    static inline m128d ceilsd(m128d a, m128d b) {
        return _mm_ceil_sd(a, b);
    }
    static inline m128d ceilpd(m128d a) {
        return _mm_ceil_pd(a);
    }
    static inline m128 floorss(m128 a, m128 b) {
        return _mm_floor_ss(a, b);
    }
    static inline m128 floorps(m128 a) {
        return _mm_floor_ps(a);
    }
    static inline m128d floorsd(m128d a, m128d b) {
        return _mm_floor_sd(a, b);
    }
    static inline m128d floorpd(m128d a) {
        return _mm_floor_pd(a);
    }
    template<int N>
    static inline m128 roundss(m128 a, m128 b) {
        return _mm_round_ss(a, b, N);
    }
    template<int N>
    static inline m128 roundps(m128 a) {
        return _mm_round_ps(a, N);
    }
    template<int N>
    static inline m128d roundsd(m128d a, m128d b) {
        return _mm_round_sd(a, b, N);
    }
    template<int N>
    static inline m128d roundpd(m128d a) {
        return _mm_round_pd(a, N);
    }

    // comparison
    static inline m128i pcmpeqq(m128i a, m128i b) {
        return _mm_cmpeq_epi64(a, b);
    }
    static inline int ptestc(m128i a, m128i b) {
        return _mm_testc_si128(a, b);
    }
    static inline int ptestz(m128i a, m128i b) {
        return _mm_testz_si128(a, b);
    }
    static inline int ptestnz(m128i a, m128i b) {
        return _mm_testnzc_si128(a, b);
    }

    // maximum, minimum
    static inline m128i pmaxsb(m128i a, m128i b) {
        return _mm_max_epi8(a, b);
    }
    static inline m128i pmaxsd(m128i a, m128i b) {
        return _mm_max_epi32(a, b);
    }
    static inline m128i pmaxuw(m128i a, m128i b) {
        return _mm_max_epu16(a, b);
    }
    static inline m128i pmaxud(m128i a, m128i b) {
        return _mm_max_epu32(a, b);
    }
    static inline m128i pminsb(m128i a, m128i b) {
        return _mm_min_epi8(a, b);
    }
    static inline m128i pminsd(m128i a, m128i b) {
        return _mm_min_epi32(a, b);
    }
    static inline m128i pminuw(m128i a, m128i b) {
        return _mm_min_epu16(a, b);
    }
    static inline m128i pminud(m128i a, m128i b) {
        return _mm_min_epu32(a, b);
    }
    static inline m128i phminposuw(m128i a) {
        return _mm_minpos_epu16(a);
    }

}

#endif
