#ifndef SIMD_SSE3_HPP
#define SIMD_SSE3_HPP

#include "./sse2.hpp"

#include <pmmintrin.h>

namespace simd {

    // move
    static inline m128i lddqu(const m128i * a) {
        return _mm_lddqu_si128( a);
    }
    static inline m128d movddup(const double * a) {
        return _mm_loaddup_pd(a);
    }
    static inline m128d movddup(m128d a) {
        return _mm_movedup_pd(a);
    }
    static inline m128 movshdup(m128 a) {
        return _mm_movehdup_ps(a);
    }
    static inline m128 movsldup(m128 a) {
        return _mm_moveldup_ps(a);
    }

    // arithmetic
    static inline m128 addsubps(m128 a, m128 b) {
        return _mm_addsub_ps(a, b);
    }
    static inline m128d addsubpd(m128d a, m128d b) {
        return _mm_addsub_pd(a, b);
    }
    static inline m128 haddps(m128 a, m128 b) {
        return _mm_hadd_ps(a, b);
    }
    static inline m128d haddpd(m128d a, m128d b) {
        return _mm_hadd_pd(a, b);
    }
    static inline m128 hsubps(m128 a, m128 b) {
        return _mm_hsub_ps(a, b);
    }
    static inline m128d hsubpd(m128d a, m128d b) {
        return _mm_hsub_pd(a, b);
    }

}

#endif
