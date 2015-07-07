#ifndef SIMD_MMX_HPP
#define SIMD_MMX_HPP

#include <cstdint>
#include <mmintrin.h>

namespace simd {

    using m64 = __m64;

    // zero
    static inline m64 zero_m64() {
        return _mm_setzero_si64();
    }

    // load
    static inline m64 load_m64(char a, char b, char c, char d, char e, char f, char g, char h) {
        return _mm_set_pi8(h, g, f, e, d, c, b, a);
    }
    static inline m64 load_m64(int8_t a, int8_t b, int8_t c, int8_t d, int8_t e, int8_t f, int8_t g, int8_t h) {
        return _mm_set_pi8(h, g, f, e, d, c, b, a);
    }
    static inline m64 load_m64(uint8_t a, uint8_t b, uint8_t c, uint8_t d, uint8_t e, uint8_t f, uint8_t g, uint8_t h) {
        return _mm_set_pi8(h, g, f, e, d, c, b, a);
    }
    static inline m64 load_m64(int16_t a, int16_t b, int16_t c, int16_t d) {
        return _mm_set_pi16(d, c, b, a);
    }
    static inline m64 load_m64(uint16_t a, uint16_t b, uint16_t c, uint16_t d) {
        return _mm_set_pi16(d, c, b, a);
    }
    static inline m64 load_m64(int32_t a, int32_t b) {
        return _mm_set_pi32(b, a);
    }
    static inline m64 load_m64(uint32_t a, uint32_t b) {
        return _mm_set_pi32(b, a);
    }
    static inline m64 load_m64(int32_t a) {
        return _mm_cvtsi32_si64(a);
    }
    static inline m64 load_m64(uint32_t a) {
        return _mm_cvtsi32_si64(a);
    }
    static inline m64 load_m64(int64_t a) {
        return _mm_cvtsi64_m64(a);
    }
    static inline m64 load_m64(uint64_t a) {
        return _mm_cvtsi64_m64(a);
    }


    // broadcast
    static inline m64 broadcast_m64(char a) {
        return _mm_set1_pi8(a);
    }
    static inline m64 broadcast_m64(int8_t a) {
        return _mm_set1_pi8(a);
    }
    static inline m64 broadcast_m64(uint8_t a) {
        return _mm_set1_pi8(a);
    }
    static inline m64 broadcast_m64(int16_t a) {
        return _mm_set1_pi16(a);
    }
    static inline m64 broadcast_m64(uint16_t a) {
        return _mm_set1_pi16(a);
    }
    static inline m64 broadcast_m64(int32_t a) {
        return _mm_set1_pi32(a);
    }
    static inline m64 broadcast_m64(uint32_t a) {
        return _mm_set1_pi32(a);
    }

    // store
    static inline int32_t m64_to_i32(m64 a) {
        return _mm_cvtsi64_si32(a);
    }
    static inline int64_t m64_to_i64(m64 a) {
        return _mm_cvtm64_si64(a);
    }

    // exit multimedia state
    static inline void emms() {
        _mm_empty();
    }

    // pack
    static inline m64 packsswb(m64 a, m64 b) {
        return _mm_packs_pi16(a, b);
    }
    static inline m64 packssdw(m64 a, m64 b) {
        return _mm_packs_pi32(a, b);
    }
    static inline m64 packuswb(m64 a, m64 b) {
        return _mm_packs_pu16(a, b);
    }

    // unpack
    static inline m64 punpckhbw(m64 a, m64 b) {
        return _mm_unpackhi_pi8(a, b);
    }
    static inline m64 punpckhwd(m64 a, m64 b) {
        return _mm_unpackhi_pi16(a, b);
    }
    static inline m64 punpckhdq(m64 a, m64 b) {
        return _mm_unpackhi_pi32(a, b);
    }
    static inline m64 punpcklbw(m64 a, m64 b) {
        return _mm_unpacklo_pi8(a, b);
    }
    static inline m64 punpcklwd(m64 a, m64 b) {
        return _mm_unpacklo_pi16(a, b);
    }
    static inline m64 punpckldq(m64 a, m64 b) {
        return _mm_unpacklo_pi32(a, b);
    }

    // addition
    static inline m64 paddb(m64 a, m64 b) {
        return _mm_add_pi8(a, b);
    }
    static inline m64 paddd(m64 a, m64 b) {
        return _mm_add_pi32(a, b);
    }
    static inline m64 paddw(m64 a, m64 b) {
        return _mm_add_pi16(a, b);
    }
    static inline m64 paddsb(m64 a, m64 b) {
        return _mm_adds_pi8(a, b);
    }
    static inline m64 paddsw(m64 a, m64 b) {
        return _mm_adds_pi16(a, b);
    }
    static inline m64 paddusb(m64 a, m64 b) {
        return _mm_adds_pu8(a, b);
    }
    static inline m64 paddusw(m64 a, m64 b) {
        return _mm_adds_pu16(a, b);
    }

    // subtraction
    static inline m64 psubb(m64 a, m64 b) {
        return _mm_sub_pi8(a, b);
    }
    static inline m64 psubw(m64 a, m64 b) {
        return _mm_sub_pi16(a, b);
    }
    static inline m64 psubd(m64 a, m64 b) {
        return _mm_sub_pi32(a, b);
    }
    static inline m64 psubsb(m64 a, m64 b) {
        return _mm_subs_pi8(a, b);
    }
    static inline m64 psubsw(m64 a, m64 b) {
        return _mm_subs_pi16(a, b);
    }
    static inline m64 psubusb(m64 a, m64 b) {
        return _mm_subs_pu8(a, b);
    }
    static inline m64 psubusw(m64 a, m64 b) {
        return _mm_subs_pu16(a, b);
    }

    // multiplication
    static inline m64 pmulhw(m64 a, m64 b) {
        return _mm_mulhi_pi16(a, b);
    }
    static inline m64 pmullw(m64 a, m64 b) {
        return _mm_mullo_pi16(a, b);
    }
    static inline m64 pmaddwd(m64 a, m64 b) {
        return _mm_madd_pi16(a, b);
    }

    // logic
    static inline m64 pand(m64 a, m64 b) {
        return _mm_and_si64(a, b);
    }
    static inline m64 pandn(m64 a, m64 b) {
        return _mm_andnot_si64(a, b);
    }
    static inline m64 por(m64 a, m64 b) {
        return _mm_or_si64(a, b);
    }
    static inline m64 pxor(m64 a, m64 b) {
        return _mm_xor_si64(a, b);
    }

    // bit shift left logical
    static inline m64 psllw(m64 a, m64 b) {
        return _mm_sll_pi16(a, b);
    }
    template<int N>
    static inline m64 psllw(m64 a) {
        return _mm_slli_pi16(a, N);
    }
    static inline m64 pslld(m64 a, m64 b) {
        return _mm_sll_pi32(a, b);
    }
    template<int N>
    static inline m64 pslld(m64 a) {
        return _mm_slli_pi32(a, N);
    }
    static inline m64 psllq(m64 a, m64 b) {
        return _mm_sll_si64(a, b);
    }
    template<int N>
    static inline m64 psllq(m64 a) {
        return _mm_slli_si64(a, N);
    }

    // bit shift right logical
    static inline m64 psrlw(m64 a, m64 b) {
        return _mm_srl_pi16(a, b);
    }
    template<int N>
    static inline m64 psrlw(m64 a) {
        return _mm_srli_pi16(a, N);
    }
    static inline m64 psrld(m64 a, m64 b) {
        return _mm_srl_pi32(a, b);
    }
    template<int N>
    static inline m64 psrld(m64 a) {
        return _mm_srli_pi32(a, N);
    }
    static inline m64 psrlq(m64 a, m64 b) {
        return _mm_srl_si64(a, b);
    }
    template<int N>
    static inline m64 psrlq(m64 a) {
        return _mm_srli_si64(a, N);
    }

    // bit shift right arithmetic
    static inline m64 psraw(m64 a, m64 b) {
        return _mm_sra_pi16(a, b);
    }
    template<int N>
    static inline m64 psraw(m64 a) {
        return _mm_srai_pi16(a, N);
    }
    static inline m64 psrad(m64 a, m64 b) {
        return _mm_sra_pi32(a, b);
    }
    template<int N>
    static inline m64 psrad(m64 a) {
        return _mm_srai_pi32(a, N);
    }

    // mask comparison
    static inline m64 pcmpeqb(m64 a, m64 b) {
        return _mm_cmpeq_pi8(a, b);
    }
    static inline m64 pcmpeqw(m64 a, m64 b) {
        return _mm_cmpeq_pi16(a, b);
    }
    static inline m64 pcmpeqd(m64 a, m64 b) {
        return _mm_cmpeq_pi32(a, b);
    }
    static inline m64 pcmpgtb(m64 a, m64 b) {
        return _mm_cmpgt_pi8(a, b);
    }
    static inline m64 pcmpgtw(m64 a, m64 b) {
        return _mm_cmpgt_pi16(a, b);
    }
    static inline m64 pcmpgtd(m64 a, m64 b) {
        return _mm_cmpgt_pi32(a, b);
    }

}

#endif
