#ifndef SIMD_SSE4_2_HPP
#define SIMD_SSE4_2_HPP

#include "sse4_1.hpp"

#include "nmmintrin.h"

namespace simd {

    // pcmpestri, pcmpestrm
    template<int i>
    static inline int pcmpestra(m128i a, int la, m128i b, int lb) {
        return _mm_cmpestra(a, la, b, lb, i);
    }
    template<int i>
    static inline int pcmpestrc(m128i a, int la, m128i b, int lb) {
        return _mm_cmpestrc(a, la, b, lb, i);
    }
    template<int i>
    static inline int pcmpestri(m128i a, int la, m128i b, int lb) {
        return _mm_cmpestri(a, la, b, lb, i);
    }
    template<int i>
    static inline int pcmpestro(m128i a, int la, m128i b, int lb) {
        return _mm_cmpestro(a, la, b, lb, i);
    }
    template<int i>
    static inline int pcmpestrs(m128i a, int la, m128i b, int lb) {
        return _mm_cmpestrs(a, la, b, lb, i);
    }
    template<int i>
    static inline int pcmpestrz(m128i a, int la, m128i b, int lb) {
        return _mm_cmpestrz(a, la, b, lb, i);
    }
    template<int i>
    static inline m128i pcmpestrm(m128i a, int la, m128i b, int lb) {
        return _mm_cmpestrm(a, la, b, lb, i);
    }

    // pcmpistri, pcmpistrm
    template<int i>
    static inline int pcmpistra(m128i a, m128i b) {
        return _mm_cmpistra(a, b, i);
    }
    template<int i>
    static inline int pcmpistrc(m128i a, m128i b) {
        return _mm_cmpistrc(a, b, i);
    }
    template<int i>
    static inline int pcmpistri(m128i a, m128i b) {
        return _mm_cmpistri(a, b, i);
    }
    template<int i>
    static inline int pcmpistro(m128i a, m128i b) {
        return _mm_cmpistro(a, b, i);
    }
    template<int i>
    static inline int pcmpistrs(m128i a, m128i b) {
        return _mm_cmpistrs(a, b, i);
    }
    template<int i>
    static inline int pcmpistrz(m128i a, m128i b) {
        return _mm_cmpistrz(a, b, i);
    }
    template<int i>
    static inline m128i pcmpistrm(m128i a, m128i b) {
        return _mm_cmpistrm(a, b, i);
    }

    // pcmpgtq
    static inline m128i pcmpgtq(m128i a, m128i b) {
        return _mm_cmpgt_epi64(a, b);
    }

    // crc32
    static inline unsigned int crc32_u8(unsigned int crc, unsigned char v) {
        return _mm_crc32_u8(crc, v);
    }
    static inline unsigned int crc32_u16(unsigned int crc, unsigned short v) {
        return _mm_crc32_u16(crc, v);
    }
    static inline unsigned int crc32_u32(unsigned int crc, unsigned int v) {
        return _mm_crc32_u32(crc, v);
    }
    static inline unsigned long long crc32_u64(unsigned long long crc, unsigned long long v) {
        return _mm_crc32_u64(crc, v);
    }

}

#endif
