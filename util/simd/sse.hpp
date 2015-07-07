#ifndef SIMD_SSE_HPP
#define SIMD_SSE_HPP

#include "mmx.hpp"

#include <xmmintrin.h>

namespace simd {

    using m128 = __m128;

    // zero
    static inline m128 zero_m128() {
        return _mm_setzero_ps();
    }

    // load
    static inline m128 load_m128_aligned(const float * a) {
        return _mm_load_ps(a);
    }
    static inline m128 load_m128(const float * a) {
        return _mm_loadu_ps(a);
    }
    static inline m128 load_m128(float a, float b, float c, float d) {
        return _mm_set_ps(d, c, b, a);
    }
    static inline m128 load_m128(float a) {
        return _mm_set_ss(a);
    }
    static inline m128 loadh_m128(m128 a, const float * b) {
        return _mm_loadh_pi(a, (m64 *) b);
    }
    static inline m128 loadl_m128(m128 a, const float * b) {
        return _mm_loadl_pi(a, (m64 *) b);
    }
    template<int N>
    static inline m64 pinsrw(m64 a, int i) {
        return _mm_insert_pi16(a, i, N);
    }

    // broadcast
    static inline m128 broadcast_m128(float a) {
        return _mm_set1_ps(a);
    }

    // store
    static inline void store_m128_aligned(float * a, m128 b) {
        _mm_store_ps(a, b);
    }
    static inline void store_m128(float * a, m128 b) {
        _mm_storeu_ps(a, b);
    }
    static inline void storeh_m128(float * a, m128 b) {
        _mm_storeh_pi((m64 *) a, b);
    }
    static inline void storel_m128(float * a, m128 b) {
        _mm_storel_pi((m64 *) a, b);
    }
    static inline float m128_to_f32(m128 a) {
        return _mm_cvtss_f32(a);
    }
    static inline void movntps(float * a, m128 b) {
        _mm_stream_ps(a, b);
    }
    static inline void movntq(m64 * a, m64 b) {
        _mm_stream_pi(a, b);
    }
    static inline void maskmovq(m64 a, m64 b, void * addr) {
        _mm_maskmove_si64(a, b, (char *) addr);
    }
    template<int N>
    static inline int32_t pextrw(m64 a) {
        return _mm_extract_pi16(a, N);
    }

    // move and shuffle
    static inline m128 movss(m128 a, m128 b) {
        return _mm_move_ss(a, b);
    }
    static inline m128 movhlps(m128 a, m128 b) {
        return _mm_movehl_ps(a, b);
    }
    static inline m128 movlhps(m128 a, m128 b) {
        return _mm_movelh_ps(a, b);
    }
    template<int i>
    static inline m128 shufps(m128 a, m128 b) {
        return _mm_shuffle_ps(a, b, i);
    }
    template<int i0, int i1, int i2, int i3>
    static inline m128 shufps(m128 a, m128 b) {
        return _mm_shuffle_ps(a, b, i0 | (i1 << 2) | (i2 << 4) | (i3 << 6));
    }
    template<int i>
    static inline m64 pshufw(m64 a) {
        return _mm_shuffle_pi16(a, i);
    }
    template<int i0, int i1, int i2, int i3>
    static inline m64 pshufw(m64 a) {
        return _mm_shuffle_pi16(a, i0 | (i1 << 2) | (i2 << 4) | (i3 << 6));
    }
    static inline int movmskps(m128 a) {
        return _mm_movemask_ps(a);
    }
    static inline int pmovmskb(m64 a) {
        return _mm_movemask_pi8(a);
    }

    // unpack
    static inline m128 unpckhps(m128 a, m128 b) {
        return _mm_unpackhi_ps(a, b);
    }
    static inline m128 unpcklps(m128 a, m128 b) {
        return _mm_unpackhi_ps(a, b);
    }

    // addition
    static inline m128 addss(m128 a, m128 b) {
        return _mm_add_ss(a, b);
    }
    static inline m128 addps(m128 a, m128 b) {
        return _mm_add_ps(a, b);
    }
    static inline m64 psadbw(m64 a, m64 b) {
        return _mm_sad_pu8(a, b);
    }

    // subtraction
    static inline m128 subss(m128 a, m128 b) {
        return _mm_sub_ss(a, b);
    }
    static inline m128 subps(m128 a, m128 b) {
        return _mm_sub_ps(a, b);
    }

    // multiplication
    static inline m128 mulss(m128 a, m128 b) {
        return _mm_mul_ss(a, b);
    }
    static inline m128 mulps(m128 a, m128 b) {
        return _mm_mul_ps(a, b);
    }
    static inline m64 pmulhuw(m64 a, m64 b) {
        return _mm_mulhi_pu16(a, b);
    }

    // division
    static inline m128 divss(m128 a, m128 b) {
        return _mm_div_ss(a, b);
    }
    static inline m128 divps(m128 a, m128 b) {
        return _mm_div_ps(a, b);
    }

    // square root
    static inline m128 sqrtss(m128 a) {
        return _mm_sqrt_ss(a);
    }
    static inline m128 sqrtps(m128 a) {
        return _mm_sqrt_ps(a);
    }

    // approximates
    static inline m128 rcpss(m128 a) {
        return _mm_rcp_ss(a);
    }
    static inline m128 rcpps(m128 a) {
        return _mm_rcp_ps(a);
    }
    static inline m128 rsqrtss(m128 a) {
        return _mm_rsqrt_ss(a);
    }
    static inline m128 rsqrtps(m128 a) {
        return _mm_rsqrt_ps(a);
    }

    // average
    static inline m64 pavgb(m64 a, m64 b) {
        return _mm_avg_pu8(a, b);
    }
    static inline m64 pavgw(m64 a, m64 b) {
        return _mm_avg_pu16(a, b);
    }

    // logic
    static inline m128 andps(m128 a, m128 b) {
        return _mm_and_ps(a, b);
    }
    static inline m128 andnps(m128 a, m128 b) {
        return _mm_andnot_ps(a, b);
    }
    static inline m128 orps(m128 a, m128 b) {
        return _mm_or_ps(a, b);
    }
    static inline m128 xorps(m128 a, m128 b) {
        return _mm_xor_ps(a, b);
    }

    // conversion
    static inline m128 cvtsi2ss(m128 a, int32_t b) {
        return _mm_cvt_si2ss(a, b);
    }
    static inline m128 cvtsi2ss(m128 a, int64_t b) {
        return _mm_cvtsi64_ss(a, b);
    }
    static inline int32_t cvtss2si(m128 a) {
        return _mm_cvtss_si32(a);
    }
    static inline int64_t cvtss2si64(m128 a) {
        return _mm_cvtss_si64(a);
    }
    static inline m128 cvtpb2ps(m64 a) {
        return _mm_cvtpi8_ps(a);
    }
    static inline m128 cvtpw2ps(m64 a) {
        return _mm_cvtpi16_ps(a);
    }
    static inline m128 cvtpi2ps(m128 a, m64 b) {
        return _mm_cvt_pi2ps(a, b);
    }
    static inline m128 cvtpi2ps(m64 a, m64 b) {
        return _mm_cvtpi32x2_ps(a, b);
    }
    static inline m64 cvtps2pb(m128 a) {
        return _mm_cvtps_pi8(a);
    }
    static inline m64 cvtps2pw(m128 a) {
        return _mm_cvtps_pi16(a);
    }
    static inline m64 cvtps2pi(m128 a) {
        return _mm_cvt_ps2pi(a);
    }
    static inline int32_t cvttss2si(m128 a) {
        return _mm_cvttss_si32(a);
    }
    static inline int64_t cvttss2si64(m128 a) {
        return _mm_cvttss_si64(a);
    }
    static inline m64 cvttps2pi(m128 a) {
        return _mm_cvttps_pi32(a);
    }

    // comparison
    static inline int comieqss(m128 a, m128 b) {
        return _mm_comieq_ss(a, b);
    }
    static inline int comigess(m128 a, m128 b) {
        return _mm_comige_ss(a, b);
    }
    static inline int comigtss(m128 a, m128 b) {
        return _mm_comigt_ss(a, b);
    }
    static inline int comiless(m128 a, m128 b) {
        return _mm_comile_ss(a, b);
    }
    static inline int comiltss(m128 a, m128 b) {
        return _mm_comilt_ss(a, b);
    }
    static inline int comineqss(m128 a, m128 b) {
        return _mm_comineq_ss(a, b);
    }
    static inline int ucomieqss(m128 a, m128 b) {
        return _mm_ucomieq_ss(a, b);
    }
    static inline int ucomigess(m128 a, m128 b) {
        return _mm_ucomige_ss(a, b);
    }
    static inline int ucomigtss(m128 a, m128 b) {
        return _mm_ucomigt_ss(a, b);
    }
    static inline int ucomiless(m128 a, m128 b) {
        return _mm_ucomile_ss(a, b);
    }
    static inline int ucomiltss(m128 a, m128 b) {
        return _mm_ucomilt_ss(a, b);
    }
    static inline int ucomineqss(m128 a, m128 b) {
        return _mm_ucomineq_ss(a, b);
    }

    // mask comparison
    static inline m128 cmpeqss(m128 a, m128 b) {
        return _mm_cmpeq_ss(a, b);
    }
    static inline m128 cmpeqps(m128 a, m128 b) {
        return _mm_cmpeq_ps(a, b);
    }
    static inline m128 cmpgess(m128 a, m128 b) {
        return _mm_cmpge_ss(a, b);
    }
    static inline m128 cmpgeps(m128 a, m128 b) {
        return _mm_cmpge_ps(a, b);
    }
    static inline m128 cmpgtss(m128 a, m128 b) {
        return _mm_cmpgt_ss(a, b);
    }
    static inline m128 cmpgtps(m128 a, m128 b) {
        return _mm_cmpgt_ps(a, b);
    }
    static inline m128 cmpless(m128 a, m128 b) {
        return _mm_cmple_ss(a, b);
    }
    static inline m128 cmpleps(m128 a, m128 b) {
        return _mm_cmple_ps(a, b);
    }
    static inline m128 cmpltss(m128 a, m128 b) {
        return _mm_cmplt_ss(a, b);
    }
    static inline m128 cmpltps(m128 a, m128 b) {
        return _mm_cmplt_ps(a, b);
    }
    static inline m128 cmpneqss(m128 a, m128 b) {
        return _mm_cmpneq_ss(a, b);
    }
    static inline m128 cmpneqps(m128 a, m128 b) {
        return _mm_cmpneq_ps(a, b);
    }
    static inline m128 cmpngess(m128 a, m128 b) {
        return _mm_cmpnge_ss(a, b);
    }
    static inline m128 cmpngeps(m128 a, m128 b) {
        return _mm_cmpnge_ps(a, b);
    }
    static inline m128 cmpngtss(m128 a, m128 b) {
        return _mm_cmpngt_ss(a, b);
    }
    static inline m128 cmpngtps(m128 a, m128 b) {
        return _mm_cmpngt_ps(a, b);
    }
    static inline m128 cmpnless(m128 a, m128 b) {
        return _mm_cmpnle_ss(a, b);
    }
    static inline m128 cmpnleps(m128 a, m128 b) {
        return _mm_cmpnle_ps(a, b);
    }
    static inline m128 cmpnltss(m128 a, m128 b) {
        return _mm_cmpnlt_ss(a, b);
    }
    static inline m128 cmpnltps(m128 a, m128 b) {
        return _mm_cmpnlt_ps(a, b);
    }
    static inline m128 cmpordss(m128 a, m128 b) {
        return _mm_cmpord_ss(a, b);
    }
    static inline m128 cmpordps(m128 a, m128 b) {
        return _mm_cmpord_ps(a, b);
    }
    static inline m128 cmpunordss(m128 a, m128 b) {
        return _mm_cmpunord_ss(a, b);
    }
    static inline m128 cmpunordps(m128 a, m128 b) {
        return _mm_cmpunord_ps(a, b);
    }

    // maximum, minimum
    static inline m128 maxss(m128 a, m128 b) {
        return _mm_max_ss(a, b);
    }
    static inline m128 maxps(m128 a, m128 b) {
        return _mm_max_ps(a, b);
    }
    static inline m128 minss(m128 a, m128 b) {
        return _mm_min_ss(a, b);
    }
    static inline m128 minps(m128 a, m128 b) {
        return _mm_min_ps(a, b);
    }
    static inline m64 pmaxub(m64 a, m64 b) {
        return _mm_max_pu8(a, b);
    }
    static inline m64 pmaxsw(m64 a, m64 b) {
        return _mm_max_pi16(a, b);
    }
    static inline m64 pminub(m64 a, m64 b) {
        return _mm_min_pu8(a, b);
    }
    static inline m64 pminsw(m64 a, m64 b) {
        return _mm_min_pi16(a, b);
    }

    // other
    static inline uint32_t stmxcsr() {
        return _mm_getcsr();
    }
    static inline void ldmxcsr(uint32_t a) {
        return _mm_setcsr(a);
    }
    static inline void prefetcht0(const char * p) {
        _mm_prefetch(p, _MM_HINT_T0);
    }
    static inline void prefetcht1(const char * p) {
        _mm_prefetch(p, _MM_HINT_T1);
    }
    static inline void prefetcht2(const char * p) {
        _mm_prefetch(p, _MM_HINT_T2);
    }
    static inline void prefetchnta(const char * p) {
        _mm_prefetch(p, _MM_HINT_NTA);
    }
    static inline void sfence() {
        _mm_sfence();
    }

    // constants
    /*static const m128 zerosf   = broadcast_m128(0.0f);
    static const m128 halfsf   = broadcast_m128(0.5f);
    static const m128 onesf    = broadcast_m128(1.0f);
    static const m128 twosf    = broadcast_m128(2.0f);
    static const m128 negonesf = broadcast_m128(-1.0f);*/

}

#endif
