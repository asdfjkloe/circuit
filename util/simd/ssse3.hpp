#ifndef SIMD_SSSE3_HPP
#define SIMD_SSSE3_HPP

#include "./sse3.hpp"

#include <tmmintrin.h>

namespace simd {

    // byte shuffle
    static inline m64 pshufb(m64 a, m64 b) {
        return _mm_shuffle_pi8(a, b);
    }
    static inline m128i pshufb(m128i a, m128i b) {
        return _mm_shuffle_epi8(a, b);
    }

    // align
    template<int N>
    static inline m64 palignr(m64 a, m64 b) {
        return _mm_alignr_pi8(a, b, N);
    }
    template<int N>
    static inline m128i palignr(m128i a, m128i b) {
        return _mm_alignr_epi8(a, b, N);
    }

    // absolute value
    static inline m64 pabsb(m64 a) {
        return _mm_abs_pi8(a);
    }
    static inline m128i pabsb(m128i a) {
        return _mm_abs_epi8(a);
    }
    static inline m64 pabsw(m64 a) {
        return _mm_abs_pi16(a);
    }
    static inline m128i pabsw(m128i a) {
        return _mm_abs_epi16(a);
    }
    static inline m64 pabsd(m64 a) {
        return _mm_abs_pi32(a);
    }
    static inline m128i pabsd(m128i a) {
        return _mm_abs_epi32(a);
    }

    // arithmetic
    static inline m64 phaddd(m64 a, m64 b) {
        return _mm_hadd_pi32(a, b);
    }
    static inline m128i phaddd(m128i a, m128i b) {
        return _mm_hadd_epi32(a, b);
    }
    static inline m64 phaddsw(m64 a, m64 b) {
        return _mm_hadds_pi16(a, b);
    }
    static inline m128i phaddsw(m128i a, m128i b) {
        return _mm_hadds_epi16(a, b);
    }
    static inline m64 phaddw(m64 a, m64 b) {
        return _mm_hadd_pi16(a, b);
    }
    static inline m128i phaddw(m128i a, m128i b) {
        return _mm_hadd_epi16(a, b);
    }
    static inline m64 phsubd(m64 a, m64 b) {
        return _mm_hsub_pi32(a, b);
    }
    static inline m128i phsubd(m128i a, m128i b) {
        return _mm_hsub_epi32(a, b);
    }
    static inline m64 phsubsw(m64 a, m64 b) {
        return _mm_hsubs_pi16(a, b);
    }
    static inline m128i phsubsw(m128i a, m128i b) {
        return _mm_hsubs_epi16(a, b);
    }
    static inline m64 phsubw(m64 a, m64 b) {
        return _mm_hsub_pi16(a, b);
    }
    static inline m128i phsubw(m128i a, m128i b) {
        return _mm_hsub_epi16(a, b);
    }
    static inline m64 pmaddubsw(m64 a, m64 b) {
        return _mm_maddubs_pi16(a, b);
    }
    static inline m128i pmaddubsw(m128i a, m128i b) {
        return _mm_maddubs_epi16(a, b);
    }
    static inline m64 pmulhrsw(m64 a, m64 b) {
        return _mm_mulhrs_pi16(a, b);
    }
    static inline m128i pmulhrsw(m128i a, m128i b) {
        return _mm_mulhrs_epi16(a, b);
    }

    // change sign
    static inline m64 psignb(m64 a, m64 b) {
        return _mm_sign_pi8(a, b);
    }
    static inline m128i psignb(m128i a, m128i b) {
        return _mm_sign_epi8(a, b);
    }
    static inline m64 psignd(m64 a, m64 b) {
        return _mm_sign_pi32(a, b);
    }
    static inline m128i psignd(m128i a, m128i b) {
        return _mm_sign_epi32(a, b);
    }
    static inline m64 psignw(m64 a, m64 b) {
        return _mm_sign_pi16(a, b);
    }
    static inline m128i psignw(m128i a, m128i b) {
        return _mm_sign_epi16(a, b);
    }

}

#endif
