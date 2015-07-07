#ifndef SIMD_SSE2_HPP
#define SIMD_SSE2_HPP

#include "./sse.hpp"

#include <emmintrin.h>

namespace simd {

    using m128d = __m128d;
    using m128i = __m128i;

    // cast
    static inline m128d m128_to_m128d(m128 a) {
        return (m128d) a;
    }
    static inline m128i m128_to_m128i(m128 a) {
        return (m128i) a;
    }
    static inline m128 m128d_to_m128(m128d a) {
        return (m128) a;
    }
    static inline m128i m128d_to_m128i(m128d a) {
        return (m128i) a;
    }
    static inline m128 m128i_to_m128(m128i a) {
        return (m128) a;
    }
    static inline m128d m128i_to_m128d(m128i a) {
        return (m128d) a;
    }

    // zero
    static inline m128d zero_m128d() {
        return _mm_setzero_pd();
    }
    static inline m128i zero_m128i() {
        return _mm_setzero_si128();
    }

    //load
    static inline m128d load_m128d_aligned(const double * a) {
        return _mm_load_pd(a);
    }
    static inline m128d load_m128d(const double * a) {
        return _mm_loadu_pd(a);
    }
    static inline m128i load_m128i_aligned(const m128i * a) {
        return _mm_load_si128(a);
    }
    static inline m128i load_m128i(const m128i * a) {
        return _mm_loadu_si128(a);
    }
    static inline m128d load_m128d(double a, double b) {
        return _mm_set_pd(b, a);
    }
    static inline m128i load_m128i(char a, char b, char c, char d, char e, char f, char g, char h,
                                   char i, char j, char k, char l, char m, char n, char o, char p) {
        return _mm_set_epi8(p, o, n, m, l, k, j, i, h, g, f, e, d, c, b, a);
    }
    static inline m128i load_m128i(int16_t a, int16_t b, int16_t c, int16_t d, int16_t e, int16_t f, int16_t g, int16_t h) {
        return _mm_set_epi16(h, g, f, e, d, c, b, a);
    }
    static inline m128i load_m128i(uint16_t a, uint16_t b, uint16_t c, uint16_t d, uint16_t e, uint16_t f, uint16_t g, uint16_t h) {
        return _mm_set_epi16(h, g, f, e, d, c, b, a);
    }
    static inline m128i load_m128i(int32_t a, int32_t b, int32_t c, int32_t d) {
        return _mm_set_epi32(d, c, b, a);
    }
    static inline m128i load_m128i(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
        return _mm_set_epi32(d, c, b, a);
    }
    static inline m128i load_m128i(int64_t a, int64_t b) {
        return _mm_set_epi64x(b, a);
    }
    static inline m128i load_m128i(uint64_t a, uint64_t b) {
        return _mm_set_epi64x(b, a);
    }
    static inline m128i load_m128i(m64 a, m64 b) {
        return _mm_set_epi64(b, a);
    }
    static inline m128d load_m128d(double a) {
        return _mm_set_sd(a);
    }
    static inline m128i load_m128i(int32_t a) {
        return _mm_cvtsi32_si128(a);
    }
    static inline m128i load_m128i(uint32_t a) {
        return _mm_cvtsi32_si128(a);
    }
    static inline m128i load_m128i(int64_t a) {
        return _mm_cvtsi64_si128(a);
    }
    static inline m128i load_m128i(uint64_t a) {
        return _mm_cvtsi64_si128(a);
    }
    static inline m128d loadh_m128d(m128d a, const double * b) {
        return _mm_loadh_pd((__m128d) a, b);
    }
    static inline m128d loadl_m128d(m128d a, const double * b) {
        return _mm_loadl_pd((__m128d) a, b);
    }
    template<int i>
    static inline m128i pinsrw(m128i a, int16_t b) {
        return _mm_insert_epi16((__m128i) a, b, i);
    }

    //broadcast
    static inline m128i broadcast_m128i(char a) {
        return _mm_set1_epi8(a);
    }
    static inline m128i broadcast_m128i(int8_t a) {
        return _mm_set1_epi8(a);
    }
    static inline m128i broadcast_m128i(uint8_t a) {
        return _mm_set1_epi8(a);
    }
    static inline m128i broadcast_m128i(int16_t a) {
        return _mm_set1_epi16(a);
    }
    static inline m128i broadcast_m128i(uint16_t a) {
        return _mm_set1_epi16(a);
    }
    static inline m128i broadcast_m128i(int32_t a) {
        return _mm_set1_epi32(a);
    }
    static inline m128i broadcast_m128i(uint32_t a) {
        return _mm_set1_epi32(a);
    }
    static inline m128i broadcast_m128i(int64_t a) {
        return _mm_set1_epi64x(a);
    }
    static inline m128i broadcast_m128i(uint64_t a) {
        return _mm_set1_epi64x(a);
    }
    static inline m128i broadcast_m128i(m64 a) {
        return _mm_set1_epi64(a);
    }
    static inline m128d broadcast_m128d(double a) {
        return _mm_set1_pd(a);
    }

    //store
    static inline void store_m128d_aligned(double * a, m128d b) {
        _mm_store_pd(a, b);
    }
    static inline void store_m128d(double * a, m128d b) {
        _mm_storeu_pd(a, b);
    }
    static inline void store_m128i_aligned(m128i * a, m128i b) {
        _mm_store_si128(a, b);
    }
    static inline void store_m128i(m128i * a, m128i b) {
        _mm_storeu_si128(a, b);
    }
    static inline void storeh_m128d(double * a, m128d b) {
        _mm_storeh_pd(a, b);
    }
    static inline void storel_m128d(double * a, m128d b) {
        _mm_storel_pd(a, b);
    }
    static inline double m128d_to_f64(m128d a) {
        return _mm_cvtsd_f64(a);
    }
    static inline int32_t m128i_to_i32(m128i a) {
        return _mm_cvtsi128_si32(a);
    }
    static inline int64_t m128i_to_i64(m128i a) {
        return _mm_cvtsi128_si64(a);
    }
    static inline void movntpd(double * a, m128d b) {
        _mm_stream_pd(a, b);
    }
    static inline void movntdq(m128i * a, m128i b) {
        _mm_stream_si128(a, b);
    }
    static inline void movnti(int32_t * a, int32_t b) {
        _mm_stream_si32(a, b);
    }
    static inline void movnti(int64_t * a, int64_t b) {
        _mm_stream_si64((long long *)a, b);
    }
    static inline void maskmovdqu(m128i a, m128i b, char * addr) {
        return _mm_maskmoveu_si128(a, b, addr);
    }
    template<int i>
    static inline int16_t pextrw(m128i a) {
        return _mm_extract_epi16(a, i);
    }

    //move and shuffle
    static inline m128d movsd(m128d a, m128d b) {
        return _mm_move_sd(a, b);
    }
    static inline m128i movq(m128i a) {
        return _mm_move_epi64(a);
    }
    static inline m64 movdq2q(m128i a) {
        return _mm_movepi64_pi64(a);
    }
    static inline m128i movq2dq(m64 a) {
        return _mm_movpi64_epi64(a);
    }
    template<int i>
    static inline m128d shufpd(m128d a, m128d b) {
        return _mm_shuffle_pd(a, b, i);
    }
    template<int i0, int i1>
    static inline m128d shufpd(m128d a, m128d b) {
        return _mm_shuffle_pd(a, b, i0 | (i1 << 1));
    }
    template<int i>
    static inline m128i pshufd(m128i a) {
        return _mm_shuffle_epi32(a, i);
    }
    template<int i0, int i1, int i2, int i3>
    static inline m128i pshufd(m128i a) {
        return _mm_shuffle_epi32(a, i0 | (i1 << 2) | (i2 << 4) | (i3 << 6));
    }
    template<int i>
    static inline m128i pshufhw(m128i a) {
        return _mm_shufflehi_epi16(a, i);
    }
    template<int i0, int i1, int i2, int i3>
    static inline m128i pshufhw(m128i a) {
        return _mm_shufflehi_epi16(a, i0 | (i1 << 2) | (i2 << 4) | (i3 << 6));
    }
    template<int i>
    static inline m128i pshuflw(m128i a) {
        return _mm_shufflelo_epi16(a, i);
    }
    template<int i0, int i1, int i2, int i3>
    static inline m128i pshuflw(m128i a) {
        return _mm_shufflelo_epi16(a, i0 | (i1 << 2) | (i2 << 4) | (i3 << 6));
    }
    static inline int movmskpd(m128d a) {
        return _mm_movemask_pd(a);
    }
    static inline int pmovmskb(m128i a) {
        return _mm_movemask_epi8(a);
    }

    //pack
    static inline m128i packsswb(m128i a, m128i b) {
        return _mm_packs_epi16(a, b);
    }
    static inline m128i packssdw(m128i a, m128i b) {
        return _mm_packs_epi32(a, b);
    }
    static inline m128i packuswb(m128i a, m128i b) {
        return _mm_packus_epi16(a, b);
    }

    //unpack
    static inline m128d unpckhpd(m128d a, m128d b) {
        return _mm_unpackhi_pd(a, b);
    }
    static inline m128d unpcklpd(m128d a, m128d b) {
        return _mm_unpackhi_pd(a, b);
    }
    static inline m128i punpckhbw(m128i a, m128i b) {
        return _mm_unpackhi_epi8(a, b);
    }
    static inline m128i punpckhdq(m128i a, m128i b) {
        return _mm_unpackhi_epi32(a, b);
    }
    static inline m128i punpckhqdq(m128i a, m128i b) {
        return _mm_unpackhi_epi64(a, b);
    }
    static inline m128i punpckhwd(m128i a, m128i b) {
        return _mm_unpackhi_epi16(a, b);
    }
    static inline m128i punpcklbw(m128i a, m128i b) {
        return _mm_unpacklo_epi8(a, b);
    }
    static inline m128i punpckldq(m128i a, m128i b) {
        return _mm_unpacklo_epi32(a, b);
    }
    static inline m128i punpcklqdq(m128i a, m128i b) {
        return _mm_unpacklo_epi64(a, b);
    }
    static inline m128i punpcklwd(m128i a, m128i b) {
        return _mm_unpacklo_epi16(a, b);
    }

    //addition
    static inline m128d addsd(m128d a, m128d b) {
        return _mm_add_sd(a, b);
    }
    static inline m128d addpd(m128d a, m128d b) {
        return _mm_add_pd(a, b);
    }
    static inline m128i paddb(m128i a, m128i b) {
        return _mm_add_epi8(a, b);
    }
    static inline m128i paddw(m128i a, m128i b) {
        return _mm_add_epi16(a, b);
    }
    static inline m128i paddd(m128i a, m128i b) {
        return _mm_add_epi32(a, b);
    }
    static inline m64 paddq(m64 a, m64 b) {
        return _mm_add_si64(a, b);
    }
    static inline m128i paddq(m128i a, m128i b) {
        return _mm_add_epi64(a, b);
    }
    static inline m128i paddsb(m128i a, m128i b) {
        return _mm_adds_epi8(a, b);
    }
    static inline m128i paddsw(m128i a, m128i b) {
        return _mm_adds_epi16(a, b);
    }
    static inline m128i paddusb(m128i a, m128i b) {
        return _mm_adds_epu8(a, b);
    }
    static inline m128i paddusw(m128i a, m128i b) {
        return _mm_adds_epu16(a, b);
    }
    static inline m128i psadbw(m128i a, m128i b) {
        return _mm_sad_epu8(a, b);
    }

    //subtraction
    static inline m128d subsd(m128d a, m128d b) {
        return _mm_sub_sd(a, b);
    }
    static inline m128d subpd(m128d a, m128d b) {
        return _mm_sub_pd(a, b);
    }
    static inline m128i psubb(m128i a, m128i b) {
        return _mm_sub_epi8(a, b);
    }
    static inline m128i psubw(m128i a, m128i b) {
        return _mm_sub_epi16(a, b);
    }
    static inline m128i psubd(m128i a, m128i b) {
        return _mm_sub_epi32(a, b);
    }
    static inline m64 psubq(m64 a, m64 b) {
        return _mm_sub_si64(a, b);
    }
    static inline m128i psubq(m128i a, m128i b) {
        return _mm_sub_epi64(a, b);
    }
    static inline m128i psubsb(m128i a, m128i b) {
        return _mm_subs_epi8(a, b);
    }
    static inline m128i psubsw(m128i a, m128i b) {
        return _mm_subs_epi16(a, b);
    }
    static inline m128i psubusb(m128i a, m128i b) {
        return _mm_subs_epu8(a, b);
    }
    static inline m128i psubusw(m128i a, m128i b) {
        return _mm_subs_epu16(a, b);
    }

    //multiplication
    static inline m128d mulsd(m128d a, m128d b) {
        return _mm_mul_sd(a, b);
    }
    static inline m128d mulpd(m128d a, m128d b) {
        return _mm_mul_pd(a, b);
    }
    static inline m128i pmulhw(m128i a, m128i b) {
        return _mm_mulhi_epi16(a, b);
    }
    static inline m128i pmullw(m128i a, m128i b) {
        return _mm_mullo_epi16(a, b);
    }
    static inline m128i pmulhuw(m128i a, m128i b) {
        return _mm_mulhi_epu16(a, b);
    }
    static inline m64 pmuludq(m64 a, m64 b) {
        return _mm_mul_su32(a, b);
    }
    static inline m128i pmuludq(m128i a, m128i b) {
        return _mm_mul_epu32(a, b);
    }
    static inline m128i pmaddwd(m128i a, m128i b) {
        return _mm_madd_epi16(a, b);
    }

    //division
    static inline m128d divsd(m128d a, m128d b) {
        return _mm_div_sd(a, b);
    }
    static inline m128d divpd(m128d a, m128d b) {
        return _mm_div_pd(a, b);
    }

    //square root
    static inline m128d sqrtsd(m128d a, m128d b) {
        //two instructions!!
        return _mm_sqrt_sd(a, b);
    }
    static inline m128d sqrtpd(m128d a) {
        return _mm_sqrt_pd(a);
    }

    //average
    static inline m128i pavgb(m128i a, m128i b) {
        return _mm_avg_epu8(a, b);
    }
    static inline m128i pavgw(m128i a, m128i b) {
        return _mm_avg_epu16(a, b);
    }

    //logic
    static inline m128d andpd(m128d a, m128d b) {
        return _mm_and_pd(a, b);
    }
    static inline m128d andnpd(m128d a, m128d b) {
        return _mm_andnot_pd(a, b);
    }
    static inline m128d orpd(m128d a, m128d b) {
        return _mm_or_pd(a, b);
    }
    static inline m128d xorpd(m128d a, m128d b) {
        return _mm_xor_pd(a, b);
    }
    static inline m128i pand(m128i a, m128i b) {
        return _mm_and_si128(a, b);
    }
    static inline m128i pandn(m128i a, m128i b) {
        return _mm_andnot_si128(a, b);
    }
    static inline m128i por(m128i a, m128i b) {
        return _mm_or_si128(a, b);
    }
    static inline m128i pxor(m128i a, m128i b) {
        return _mm_xor_si128(a, b);
    }

    //bit shift
    template<int i>
    static inline m128i pslld(m128i a) {
        return _mm_slli_epi32(a, i);
    }
    static inline m128i pslld(m128i a, m128i b) {
        return _mm_sll_epi32(a, b);
    }
    template<int i>
    static inline m128i pslldq(m128i a) {
        return _mm_slli_si128((__m128i)a, i);
    }
    template<int i>
    static inline m128i psllq(m128i a) {
        return _mm_slli_epi64(a, i);
    }
    static inline m128i psllq(m128i a, m128i b) {
        return _mm_sll_epi64(a, b);
    }
    template<int i>
    static inline m128i psllw(m128i a) {
        return _mm_slli_epi16(a, i);
    }
    static inline m128i psllw(m128i a, m128i b) {
        return _mm_sll_epi16(a, b);
    }
    template<int i>
    static inline m128i psrad(m128i a) {
        return _mm_srai_epi32(a, i);
    }
    static inline m128i psrad(m128i a, m128i b) {
        return _mm_sra_epi32(a, b);
    }
    template<int i>
    static inline m128i psraw(m128i a) {
        return _mm_srai_epi16(a, i);
    }
    static inline m128i psraw(m128i a, m128i b) {
        return _mm_sra_epi16(a, b);
    }
    template<int i>
    static inline m128i psrld(m128i a) {
        return _mm_srli_epi32(a, i);
    }
    static inline m128i psrld(m128i a, m128i b) {
        return _mm_srl_epi32(a, b);
    }
    template<int i>
    static inline m128i psrldq(m128i a) {
        return _mm_srli_si128((__m128i)a, i);
    }
    template<int i>
    static inline m128i psrlq(m128i a) {
        return _mm_srli_epi64(a, i);
    }
    static inline m128i psrlq(m128i a, m128i b) {
        return _mm_srl_epi64(a, b);
    }
    template<int i>
    static inline m128i psrlw(m128i a) {
        return _mm_srli_epi16(a, i);
    }
    static inline m128i psrlw(m128i a, m128i b) {
        return _mm_srl_epi16(a, b);
    }

    //conversion
    static inline m128d cvtdq2pd(m128i a) {
        return _mm_cvtepi32_pd(a);
    }
    static inline m128 cvtdq2ps(m128i a) {
        return _mm_cvtepi32_ps(a);
    }
    static inline m128i cvtpd2dq(m128d a) {
        return _mm_cvtpd_epi32(a);
    }
    static inline m64 cvtpd2pi(m128d a) {
        return _mm_cvtpd_pi32(a);
    }
    static inline m128 cvtpd2ps(m128d a) {
        return _mm_cvtpd_ps(a);
    }
    static inline m128d cvtpi2pd(m64 a) {
        return _mm_cvtpi32_pd(a);
    }
    static inline m128i cvtps2dq(m128 a) {
        return _mm_cvtps_epi32(a);
    }
    static inline m128d cvtps2pd(m128 a) {
        return _mm_cvtps_pd(a);
    }
    static inline int32_t cvtsd2si(m128d a) {
        return _mm_cvtsd_si32(a);
    }
    static inline int64_t cvtsd2si64(m128d a) {
        return _mm_cvtsd_si64(a);
    }
    static inline m128 cvtsd2ss(m128 a, m128d b) {
        return _mm_cvtsd_ss(a, b);
    }
    static inline m128d cvtsi2sd(m128d a, int32_t b) {
        return _mm_cvtsi32_sd(a, b);
    }
    static inline m128d cvtsi2sd(m128d a, int64_t b) {
        return _mm_cvtsi64_sd(a, b);
    }
    static inline m128d cvtss2sd(m128d a, m128 b) {
        return _mm_cvtss_sd(a, b);
    }
    static inline m128i cvttpd2dq(m128d a) {
        return _mm_cvttpd_epi32(a);
    }
    static inline m64 cvttpd2pi(m128d a) {
        return _mm_cvttpd_pi32(a);
    }
    static inline m128i cvttps2dq(m128 a) {
        return _mm_cvttps_epi32(a);
    }
    static inline int32_t cvttsd2si(m128d a) {
        return _mm_cvttsd_si32(a);
    }
    static inline int64_t cvttsd2si64(m128d a) {
        return _mm_cvttsd_si64(a);
    }

    //comparison
    static inline int comieqsd(m128d a, m128d b) {
        return _mm_comieq_sd(a, b);
    }
    static inline int comigesd(m128d a, m128d b) {
        return _mm_comige_sd(a, b);
    }
    static inline int comigtsd(m128d a, m128d b) {
        return _mm_comigt_sd(a, b);
    }
    static inline int comilesd(m128d a, m128d b) {
        return _mm_comile_sd(a, b);
    }
    static inline int comiltsd(m128d a, m128d b) {
        return _mm_comilt_sd(a, b);
    }
    static inline int comineqsd(m128d a, m128d b) {
        return _mm_comineq_sd(a, b);
    }
    static inline int ucomieqsd(m128d a, m128d b) {
        return _mm_ucomieq_sd(a, b);
    }
    static inline int ucomigesd(m128d a, m128d b) {
        return _mm_ucomige_sd(a, b);
    }
    static inline int ucomigtsd(m128d a, m128d b) {
        return _mm_ucomigt_sd(a, b);
    }
    static inline int ucomilesd(m128d a, m128d b) {
        return _mm_ucomile_sd(a, b);
    }
    static inline int ucomiltsd(m128d a, m128d b) {
        return _mm_ucomilt_sd(a, b);
    }
    static inline int ucomineqsd(m128d a, m128d b) {
        return _mm_ucomineq_sd(a, b);
    }

    //mask comparison
    static inline m128d cmpeqsd(m128d a, m128d b) {
        return _mm_cmpeq_sd(a, b);
    }
    static inline m128d cmpeqpd(m128d a, m128d b) {
        return _mm_cmpeq_pd(a, b);
    }
    static inline m128d cmpgesd(m128d a, m128d b) {
        return _mm_cmpge_sd(a, b);
    }
    static inline m128d cmpgepd(m128d a, m128d b) {
        return _mm_cmpge_pd(a, b);
    }
    static inline m128d cmpgtsd(m128d a, m128d b) {
        return _mm_cmpgt_sd(a, b);
    }
    static inline m128d cmpgtpd(m128d a, m128d b) {
        return _mm_cmpgt_pd(a, b);
    }
    static inline m128d cmplesd(m128d a, m128d b) {
        return _mm_cmple_sd(a, b);
    }
    static inline m128d cmplepd(m128d a, m128d b) {
        return _mm_cmple_pd(a, b);
    }
    static inline m128d cmpltsd(m128d a, m128d b) {
        return _mm_cmplt_sd(a, b);
    }
    static inline m128d cmpltpd(m128d a, m128d b) {
        return _mm_cmplt_pd(a, b);
    }
    static inline m128d cmpordsd(m128d a, m128d b) {
        return _mm_cmpord_sd(a, b);
    }
    static inline m128d cmpordpd(m128d a, m128d b) {
        return _mm_cmpord_pd(a, b);
    }
    static inline m128d cmpneqsd(m128d a, m128d b) {
        return _mm_cmpneq_sd(a, b);
    }
    static inline m128d cmpneqpd(m128d a, m128d b) {
        return _mm_cmpneq_pd(a, b);
    }
    static inline m128d cmpngesd(m128d a, m128d b) {
        return _mm_cmpnge_sd(a, b);
    }
    static inline m128d cmpngepd(m128d a, m128d b) {
        return _mm_cmpnge_pd(a, b);
    }
    static inline m128d cmpngtsd(m128d a, m128d b) {
        return _mm_cmpngt_sd(a, b);
    }
    static inline m128d cmpngtpd(m128d a, m128d b) {
        return _mm_cmpngt_pd(a, b);
    }
    static inline m128d cmpnlesd(m128d a, m128d b) {
        return _mm_cmpnle_sd(a, b);
    }
    static inline m128d cmpnlepd(m128d a, m128d b) {
        return _mm_cmpnle_pd(a, b);
    }
    static inline m128d cmpnltsd(m128d a, m128d b) {
        return _mm_cmpnlt_sd(a, b);
    }
    static inline m128d cmpnltpd(m128d a, m128d b) {
        return _mm_cmpnlt_pd(a, b);
    }
    static inline m128d cmpunordsd(m128d a, m128d b) {
        return _mm_cmpunord_sd(a, b);
    }
    static inline m128d cmpunordpd(m128d a, m128d b) {
        return _mm_cmpunord_pd(a, b);
    }
    static inline m128i pcmpeqb(m128i a, m128i b) {
        return _mm_cmpeq_epi8(a, b);
    }
    static inline m128i pcmpeqw(m128i a, m128i b) {
        return _mm_cmpeq_epi16(a, b);
    }
    static inline m128i pcmpeqd(m128i a, m128i b) {
        return _mm_cmpeq_epi32(a, b);
    }
    static inline m128i pcmpgtb(m128i a, m128i b) {
        return _mm_cmpgt_epi8(a, b);
    }
    static inline m128i pcmpgtw(m128i a, m128i b) {
        return _mm_cmpgt_epi16(a, b);
    }
    static inline m128i pcmpgtd(m128i a, m128i b) {
        return _mm_cmpgt_epi32(a, b);
    }
    static inline m128i pcmpltb(m128i a, m128i b) {
        return _mm_cmplt_epi8(a, b);
    }
    static inline m128i pcmpltw(m128i a, m128i b) {
        return _mm_cmplt_epi16(a, b);
    }
    static inline m128i pcmpltd(m128i a, m128i b) {
        return _mm_cmplt_epi32(a, b);
    }

    //maximum, minimum
    static inline m128d maxsd(m128d a, m128d b) {
        return _mm_max_sd(a, b);
    }
    static inline m128d maxpd(m128d a, m128d b) {
        return _mm_max_pd(a, b);
    }
    static inline m128d minsd(m128d a, m128d b) {
        return _mm_min_sd(a, b);
    }
    static inline m128d minpd(m128d a, m128d b) {
        return _mm_min_pd(a, b);
    }
    static inline m128i pmaxub(m128i a, m128i b) {
        return _mm_max_epu8(a, b);
    }
    static inline m128i pmaxsw(m128i a, m128i b) {
        return _mm_max_epi16(a, b);
    }
    static inline m128i pminub(m128i a, m128i b) {
        return _mm_min_epu8(a, b);
    }
    static inline m128i pminsw(m128i a, m128i b) {
        return _mm_min_epi16(a, b);
    }

    //other
    static inline void clflush(const void * a) {
        _mm_clflush(a);
    }
    static inline void lfence() {
        _mm_lfence();
    }
    static inline void mfence() {
        _mm_mfence();
    }
    static inline void pause() {
        _mm_pause();
    }

    //constants
    /*static const m128d zerosd   = broadcast_m128d(0.0);
    static const m128d halfsd   = broadcast_m128d(0.5);
    static const m128d onesd    = broadcast_m128d(1.0);
    static const m128d twosd    = broadcast_m128d(2.0);
    static const m128d negonesd = broadcast_m128d(-1.0);
    static const m128i zeros    = broadcast_m128i(0x00000000);
    static const m128i ones     = broadcast_m128i(0xFFFFFFFF);
    static const m128i sgnmsk   = broadcast_m128i(0x80000000);*/

}

#endif
