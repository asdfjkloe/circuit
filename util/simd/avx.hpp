#ifndef SIMD_AVX_HPP
#define SIMD_AVX_HPP

#include "sse4_2.hpp"

#include "immintrin.h"

namespace simd {

    using m256  = __m256;
    using m256d = __m256d;
    using m256i = __m256i;

    // cast
    static inline m256 m128_to_m256(m128 a) {
        return _mm256_castps128_ps256(a);
    }
    static inline m256d m128d_to_m256d(m128d a) {
        return _mm256_castpd128_pd256(a);
    }
    static inline m256i m128i_to_m256i(m128i a) {
        return _mm256_castsi128_si256(a);
    }
    static inline m128 m256_to_m128(m256 a) {
        return _mm256_castps256_ps128(a);
    }
    static inline m256d m256_to_m256d(m256 a) {
        return _mm256_castps_pd(a);
    }
    static inline m256i m256_to_m256i(m256 a) {
        return _mm256_castps_si256(a);
    }
    static inline m128d m256d_to_m128d(m256d a) {
        return _mm256_castpd256_pd128(a);
    }
    static inline m256 m256d_to_m256(m256d a) {
        return _mm256_castpd_ps(a);
    }
    static inline m256i m256d_to_m256i(m256d a) {
        return _mm256_castpd_si256(a);
    }
    static inline m128i m256i_to_m128i(m256i a) {
        return _mm256_castsi256_si128(a);
    }
    static inline m256 m256i_to_m256(m256i a) {
        return _mm256_castsi256_ps(a);
    }
    static inline m256d m256i_to_m256d(m256i a) {
        return _mm256_castsi256_pd(a);
    }

    // zero
    static inline m256 zero_m256() {
        return _mm256_setzero_ps();
    }
    static inline m256d zero_m256d() {
        return _mm256_setzero_pd();
    }
    static inline m256i zero_m256i() {
        return _mm256_setzero_si256();
    }
    static inline void vzeroall() {
        _mm256_zeroall();
    }
    static inline void vzeroupper() {
        _mm256_zeroupper();
    }

    // load
    static inline m256 load_m256_aligned(const float * a) {
        return _mm256_load_ps(a);
    }
    static inline m256 load_m256(const float * a) {
        return _mm256_loadu_ps(a);
    }
    /*static inline m256 load_m256(const float * a, const float * b) {
        return _mm256_loadu2_m128(b, a);
    }*/
    static inline m256d load_m256d_aligned(const double * a) {
        return _mm256_load_pd(a);
    }
    static inline m256d load_m256d(const double * a) {
        return _mm256_loadu_pd(a);
    }
    /*static inline m256d load_m256d(const double * a, const double * b) {
        return _mm256_loadu2_m128d(b, a);
    }*/
    static inline m256i load_m256i_aligned(const m256i * a) {
        return _mm256_load_si256(a);
    }
    static inline m256i load_m256i(const m256i * a) {
        return _mm256_loadu_si256(a);
    }
    /*static inline m256i load_m256i(const m128i * a, const m128i * b) {
        return _mm256_loadu2_m128i(b, a);
    }*/
    static inline m256 load_m256(float a, float b, float c, float d, float e, float f, float g, float h) {
        return _mm256_set_ps(h, g, f, e, d, c, b, a);
    }
    static inline m256d load_m256d(double a, double b, double c, double d) {
        return _mm256_set_pd(d, c, b, a);
    }
    static inline m256i load_m256i(char a00, char a01, char a02, char a03, char a04, char a05, char a06, char a07,
                                   char a08, char a09, char a10, char a11, char a12, char a13, char a14, char a15,
                                   char a16, char a17, char a18, char a19, char a20, char a21, char a22, char a23,
                                   char a24, char a25, char a26, char a27, char a28, char a29, char a30, char a31) {
        return _mm256_set_epi8(a31, a30, a29, a28, a27, a26, a25, a24,
                               a23, a22, a21, a20, a19, a18, a17, a16,
                               a15, a14, a13, a12, a11, a10, a09, a08,
                               a07, a06, a05, a04, a03, a02, a01, a00);
    }
    static inline m256i load_m256i(int16_t a, int16_t b, int16_t c, int16_t d, int16_t e, int16_t f, int16_t g, int16_t h,
                                   int16_t i, int16_t j, int16_t k, int16_t l, int16_t m, int16_t n, int16_t o, int16_t p) {
        return _mm256_set_epi16(p, o, n, m, l, k, j, i, h, g, f, e, d, c, b, a);
    }
    static inline m256i load_m256i(uint16_t a, uint16_t b, uint16_t c, uint16_t d, uint16_t e, uint16_t f, uint16_t g, uint16_t h,
                                   uint16_t i, uint16_t j, uint16_t k, uint16_t l, uint16_t m, uint16_t n, uint16_t o, uint16_t p) {
        return _mm256_set_epi16(p, o, n, m, l, k, j, i, h, g, f, e, d, c, b, a);
    }
    static inline m256i load_m256i(int32_t a, int32_t b, int32_t c, int32_t d, int32_t e, int32_t f, int32_t g, int32_t h) {
        return _mm256_set_epi32(h, g, f, e, d, c, b, a);
    }
    static inline m256i load_m256i(uint32_t a, uint32_t b, uint32_t c, uint32_t d, uint32_t e, uint32_t f, uint32_t g, uint32_t h) {
        return _mm256_set_epi32(h, g, f, e, d, c, b, a);
    }
    static inline m256i load_m256i(int64_t a, int64_t b, int64_t c, int64_t d) {
        return _mm256_set_epi64x(d, c, b, a);
    }
    static inline m256i load_m256i(uint64_t a, uint64_t b, uint64_t c, uint64_t d) {
        return _mm256_set_epi64x(d, c, b, a);
    }
    static inline m128 load_m128() {
        return _mm_undefined_ps();
    }
    static inline m128d load_m128d() {
        return _mm_undefined_pd();
    }
    static inline m128i load_m128i() {
        return _mm_undefined_si128();
    }
    static inline m256 load_m256() {
        return _mm256_undefined_ps();
    }
    static inline m256d load_m256d() {
        return _mm256_undefined_pd();
    }
    static inline m256i load_m256i() {
        return _mm256_undefined_si256();
    }
    static inline m256i vlddqu(const m256i * a) {
        return _mm256_lddqu_si256(a);
    }
    static inline m128 vmaskmovps128(const float * a, m128i msk) {
        return _mm_maskload_ps(a, msk);
    }
    static inline m128d vmaskmovpd128(const double * a, m128i msk) {
        return _mm_maskload_pd(a, msk);
    }
    static inline m256 vmaskmovps(const float * a, m256i msk) {
        return _mm256_maskload_ps(a, msk);
    }
    static inline m256d vmaskmovpd(const double * a, m256i msk) {
        return _mm256_maskload_pd(a, msk);
    }

    // broadcast
    static inline m256 broadcast_m256(float a) {
        return _mm256_set1_ps(a);
    }
    static inline m256d broadcast_m256d(double a) {
        return _mm256_set1_pd(a);
    }
    static inline m256i broadcast_m256i(char a) {
        return _mm256_set1_epi8(a);
    }
    static inline m256i broadcast_m256i(int8_t a) {
        return _mm256_set1_epi8(a);
    }
    static inline m256i broadcast_m256i(uint8_t a) {
        return _mm256_set1_epi8(a);
    }
    static inline m256i broadcast_m256i(int16_t a) {
        return _mm256_set1_epi16(a);
    }
    static inline m256i broadcast_m256i(uint16_t a) {
        return _mm256_set1_epi16(a);
    }
    static inline m256i broadcast_m256i(int32_t a) {
        return _mm256_set1_epi32(a);
    }
    static inline m256i broadcast_m256i(uint32_t a) {
        return _mm256_set1_epi32(a);
    }
    static inline m256i broadcast_m256i(int64_t a) {
        return _mm256_set1_epi64x(a);
    }
    static inline m256i broadcast_m256i(uint64_t a) {
        return _mm256_set1_epi64x(a);
    }
    static inline m128 vbroadcastss128(const float * a) {
        return _mm_broadcast_ss(a);
    }
    static inline m256 vbroadcastss(const float * a) {
        return _mm256_broadcast_ss(a);
    }
    static inline m256d vbroadcastsd(const double * a) {
        return _mm256_broadcast_sd(a);
    }
    static inline m256 vbroadcastf128(const m128 * a) {
        return _mm256_broadcast_ps(a);
    }
    static inline m256d vbroadcastf128(const m128d * a) {
        return _mm256_broadcast_pd(a);
    }

    // store
    static inline void store_m256_aligned(float * a, m256 b) {
        _mm256_store_ps(a, b);
    }
    static inline void store_m256(float * a, m256 b) {
        _mm256_storeu_ps(a, b);
    }
    /*static inline void store_m256(float * a, float * b, m256 c) {
        _mm256_storeu2_m128(a, b, c);
    }*/
    static inline void store_m256d_aligned(double * a, m256d b) {
        _mm256_store_pd(a, b);
    }
    static inline void store_m256d(double * a, m256d b) {
        _mm256_storeu_pd(a, b);
    }
    /*static inline void store_m256d(double * a, double * b, m256d c) {
        _mm256_storeu2_m128d(a, b, c);
    }*/
    static inline void store_m256i_aligned(m256i * a, m256i b) {
        _mm256_store_si256(a, b);
    }
    static inline void store_m256i(m256i * a, m256i b) {
        _mm256_storeu_si256(a, b);
    }
    /*static inline void store_m256i(m128i * a, m128i * b, m256 c) {
        _mm256_storeu2_m128i(a, b, c);
    }*/
    static inline void vmaskmovps128(float * a, m128i msk, m128 b) {
        _mm_maskstore_ps(a, msk, b);
    }
    static inline void vmaskmovps(float * a, m256i msk, m256 b) {
        _mm256_maskstore_ps(a, msk, b);
    }
    static inline void vmaskmovpd128(double * a, m128i msk, m128d b) {
        _mm_maskstore_pd(a, msk, b);
    }
    static inline void vmaskmovpd(double * a, m256i msk, m256d b) {
        _mm256_maskstore_pd(a, msk, b);
    }
    static inline void vmovntps(float * a, m256 b) {
        _mm256_stream_ps(a, b);
    }
    static inline void vmovntpd(double * a, m256d b) {
        _mm256_stream_pd(a, b);
    }
    static inline void vmovntdq(m256i * a, m256i b) {
        _mm256_stream_si256(a, b);
    }

    // move and shuffle
    static inline m256 vmovsldup(m256 a) {
        return _mm256_moveldup_ps(a);
    }
    static inline m256 vmovshdup(m256 a) {
        return _mm256_movehdup_ps(a);
    }
    static inline m256d vmovddup(m256d a) {
        return _mm256_movedup_pd(a);
    }
    static inline int vmovmskps(m256 a) {
        return _mm256_movemask_ps(a);
    }
    static inline int vmovmskpd(m256d a) {
        return _mm256_movemask_pd(a);
    }
    template<int i>
    static inline m256 vshufps(m256 a, m256 b) {
        return _mm256_shuffle_ps(a, b, i);
    }
    template<int i0, int i1, int i2, int i3>
    static inline m256 vshufps(m256 a, m256 b) {
        return _mm256_shuffle_ps(a, b, i0 | (i1 << 2) | (i2 << 4) | (i3 << 6));
    }
    template<int i>
    static inline m256d vshufpd(m256d a, m256d b) {
        return _mm256_shuffle_pd(a, b, i);
    }
    template<int i0, int i1, int i2, int i3>
    static inline m256d vshufpd(m256d a, m256d b) {
        return _mm256_shuffle_pd(a, b, i0 | (i1 << 1) | (i2 << 2) | (i3 << 4));
    }
    template<int i>
    static inline m128 vpermilps(m128 a) {
        return _mm_permute_ps(a, i);
    }
    template<int i0, int i1, int i2, int i3>
    static inline m128 vpermilps(m128 a) {
        return _mm_permute_ps(a, i0 | (i1 << 2) | (i2 << 4) | (i3 << 6));
    }
    template<int i>
    static inline m256 vpermilps(m256 a) {
        return _mm256_permute_ps(a, i);
    }
    template<int i0, int i1, int i2, int i3>
    static inline m256 vpermilps(m256 a) {
        return _mm256_permute_ps(a, i0 | (i1 << 2) | (i2 << 4) | (i3 << 6));
    }
    template<int i>
    static inline m128d vpermilpd(m128d a) {
        return _mm_permute_pd(a, i);
    }
    template<int i0, int i1>
    static inline m128d vpermilpd(m128d a) {
        return _mm_permute_pd(a, i0 | (i1 << 1));
    }
    template<int i>
    static inline m256d vpermilpd(m256d a) {
        return _mm256_permute_pd(a, i);
    }
    template<int i0, int i1, int i2, int i3>
    static inline m256d vpermilpd(m256d a) {
        return _mm256_permute_pd(a, i0 | (i1 << 1) | (i2 << 2) | (i3 << 3));
    }
    template<int i>
    static inline m256 vperm2f128(m256 a, m256 b) {
        return _mm256_permute2f128_ps(a, b, i);
    }
    template<int i0, int i1>
    static inline m256 vperm2f128(m256 a, m256 b) {
        return _mm256_permute2f128_ps(a, b, i0 | (i1 << 4));
    }
    template<int i>
    static inline m256d vperm2f128(m256d a, m256d b) {
        return _mm256_permute2f128_pd(a, b, i);
    }
    template<int i0, int i1>
    static inline m256d vperm2f128(m256d a, m256d b) {
        return _mm256_permute2f128_pd(a, b, i0 | (i1 << 4));
    }
    template<int i>
    static inline m256i vperm2f128(m256i a, m256i b) {
        return _mm256_permute2f128_si256(a, b, i);
    }
    template<int i0, int i1>
    static inline m256i vperm2f128(m256i a, m256i b) {
        return _mm256_permute2f128_si256(a, b, i0 | (i1 << 4));
    }
    static inline m128 vpermilps128(m128 a, m128i b) {
        return _mm_permutevar_ps(a, b);
    }
    static inline m256 vpermilps(m256 a, m256i b) {
        return _mm256_permutevar_ps(a, b);
    }
    static inline m128d vpermilpd128(m128d a, m128i b) {
        return _mm_permutevar_pd(a, b);
    }
    static inline m256d vpermilpd(m256d a, m256i b) {
        return _mm256_permutevar_pd(a, b);
    }

    // unpack
    static inline m256 vunpcklps(m256 a, m256 b) {
        return _mm256_unpacklo_ps(a, b);
    }
    static inline m256 vunpckhps(m256 a, m256 b) {
        return _mm256_unpackhi_ps(a, b);
    }
    static inline m256d vunpcklpd(m256d a, m256d b) {
        return _mm256_unpacklo_pd(a, b);
    }
    static inline m256d vunpckhpd(m256d a, m256d b) {
        return _mm256_unpackhi_pd(a, b);
    }

    // extract and insert
    template<int i>
    static inline int8_t vpextrb(m256i a) {
        return _mm256_extract_epi8(a, i);
    }
    template<int i>
    static inline int16_t vpextrw(m256i a) {
        return _mm256_extract_epi16(a, i);
    }
    template<int i>
    static inline int32_t vpextrd(m256i a) {
        return _mm256_extract_epi32(a, i);
    }
    template<int i>
    static inline int64_t vpextrq(m256i a) {
        return _mm256_extract_epi64(a, i);
    }
    template<int i>
    static inline m128 vextractf128(m256 a) {
        return _mm256_extractf128_ps(a, i);
    }
    template<int i>
    static inline m128d vextractf128(m256d a) {
        return _mm256_extractf128_pd(a, i);
    }
    template<int i>
    static inline m128i vextractf128(m256i a) {
        return _mm256_extractf128_si256(a, i);
    }
    template<int i>
    static inline m256i vpinsrb(m256i a, int8_t b) {
        return _mm256_insert_epi8(a, b, i);
    }
    template<int i>
    static inline m256i vpinsrw(m256i a, int16_t b) {
        return _mm256_insert_epi16(a, b, i);
    }
    template<int i>
    static inline m256i vpinsrd(m256i a, int32_t b) {
        return _mm256_insert_epi32(a, b, i);
    }
    template<int i>
    static inline m256i vpinsrq(m256i a, int64_t b) {
        return _mm256_insert_epi64(a, b, i);
    }
    template<int i>
    static inline m256 vinsertf128(m256 a, m128 b) {
        return _mm256_insertf128_ps(a, b, i);
    }
    template<int i>
    static inline m256d vinsertf128(m256d a, m128d b) {
        return _mm256_insertf128_pd(a, b, i);
    }
    template<int i>
    static inline m256i vinsertf128(m256i a, m128i b) {
        return _mm256_insertf128_si256(a, b, i);
    }
    /*static inline m256 vinsertf128(m128 a, m128 b) {
        return _mm256_set_m128(b, a);
    }*/
    /*static inline m256d vinsertf128(m128d a, m128d b) {
        return _mm256_set_m128d(b, a);
    }*/
    /*static inline m256i vinsertf128(m128i a, m128i b) {
        return _mm256_set_m128i(b, a);
    }*/

    // blend
    template<int i>
    static inline m256 vblendps(m256 a, m256 b) {
        return _mm256_blend_ps (a, b, i);
    }
    template<int i>
    static inline m256d vblendpd(m256d a, m256d b) {
        return _mm256_blend_pd (a, b, i);
    }
    static inline m256 vblendvps(m256 a, m256 b, m256 msk) {
        return _mm256_blendv_ps (a, b, msk);
    }
    static inline m256d vblendvpd(m256d a, m256d b, m256d msk) {
        return _mm256_blendv_pd (a, b, msk);
    }

    // addition
    static inline m256 vaddps(m256 a, m256 b) {
        return _mm256_add_ps(a, b);
    }
    static inline m256d vaddpd(m256d a, m256d b) {
        return _mm256_add_pd(a, b);
    }
    static inline m256 vhaddps(m256 a, m256 b) {
        return _mm256_hadd_ps(a, b);
    }
    static inline m256d vhaddpd(m256d a, m256d b) {
        return _mm256_hadd_pd(a, b);
    }

    // subtraction
    static inline m256 vsubps(m256 a, m256 b) {
        return _mm256_sub_ps(a, b);
    }
    static inline m256d vsubpd(m256d a, m256d b) {
        return _mm256_sub_pd(a, b);
    }
    static inline m256 vhsubps(m256 a, m256 b) {
        return _mm256_hsub_ps(a, b);
    }
    static inline m256d vhsubpd(m256d a, m256d b) {
        return _mm256_hsub_pd(a, b);
    }

    // addsub
    static inline m256 vaddsubps(m256 a, m256 b) {
        return _mm256_addsub_ps(a, b);
    }
    static inline m256d vaddsubpd(m256d a, m256d b) {
        return _mm256_addsub_pd(a, b);
    }

    // multiplication
    static inline m256 vmulps(m256 a, m256 b) {
        return _mm256_mul_ps(a, b);
    }
    static inline m256d vmulpd(m256d a, m256d b) {
        return _mm256_mul_pd(a, b);
    }

    // division
    static inline m256 vdivps(m256 a, m256 b) {
        return _mm256_div_ps(a, b);
    }
    static inline m256d vdivpd(m256d a, m256d b) {
        return _mm256_div_pd(a, b);
    }

    // dot product
    template<int i>
    static inline m256 vdpps(m256 a, m256 b) {
        return _mm256_dp_ps(a, b, i);
    }

    // square root
    static inline m256 vsqrtps(m256 a) {
        return _mm256_sqrt_ps(a);
    }
    static inline m256d vsqrtpd(m256d a) {
        return _mm256_sqrt_pd(a);
    }

    // approximates
    static inline m256 vrcpps(m256 a) {
        return _mm256_rcp_ps(a);
    }
    static inline m256 vrsqrtps(m256 a) {
        return _mm256_rsqrt_ps(a);
    }

    // round
    static inline m256 vceilps(m256 a) {
        return _mm256_ceil_ps(a);
    }
    static inline m256d vceilpd(m256d a) {
        return _mm256_ceil_pd(a);
    }
    static inline m256 vfloorps(m256 a) {
        return _mm256_floor_ps(a);
    }
    static inline m256d vfloorpd(m256d a) {
        return _mm256_floor_pd(a);
    }
    template<int i>
    static inline m256 vroundps(m256 a) {
        return _mm256_round_ps(a, i);
    }
    template<int i>
    static inline m256d vroundpd(m256d a) {
        return _mm256_round_pd(a, i);
    }

    // logic
    static inline m256 vandps(m256 a, m256 b) {
        return _mm256_and_ps(a, b);
    }
    static inline m256d vandpd(m256d a, m256d b) {
        return _mm256_and_pd(a, b);
    }
    static inline m256 vandnps(m256 a, m256 b){
        return _mm256_andnot_ps(a, b);
    }
    static inline m256d vandnpd(m256d a, m256d b){
        return _mm256_andnot_pd(a, b);
    }
    static inline m256 vorps(m256 a, m256 b){
        return _mm256_or_ps(a, b);
    }
    static inline m256d vorpd(m256d a, m256d b){
        return _mm256_or_pd(a, b);
    }
    static inline m256 vxorps(m256 a, m256 b){
        return _mm256_xor_ps(a, b);
    }
    static inline m256d vxorpd(m256d a, m256d b){
        return _mm256_xor_pd(a, b);
    }

    // conversion
    static inline m256i vcvtps2dq(m256 a) {
        return _mm256_cvtps_epi32(a);
    }
    static inline m256d vcvtps2pd(m128 a) {
        return _mm256_cvtps_pd(a);
    }
    static inline m128 vcvtpd2ps(m256d a) {
        return _mm256_cvtpd_ps(a);
    }
    static inline m128i vcvtpd2dq(m256d a) {
        return _mm256_cvtpd_epi32(a);
    }
    static inline m256 vcvtdq2ps(m256i a) {
        return _mm256_cvtepi32_ps(a);
    }
    static inline m256d vcvtdq2pd(m128i a) {
        return _mm256_cvtepi32_pd(a);
    }
    static inline m256i vcvttps2dq(m256 a) {
        return _mm256_cvttps_epi32(a);
    }
    static inline m128i vcvttpd2dq(m256d a) {
        return _mm256_cvttpd_epi32(a);
    }


    // comparison
    template<int i>
    static inline m128 vcmpps(m128 a, m128 b) {
        return _mm_cmp_ps(a, b, i);
    }
    template<int i>
    static inline m256 vcmpps(m256 a, m256 b) {
        return _mm256_cmp_ps(a, b, i);
    }
    template<int i>
    static inline m128d vcmppd(m128d a, m128d b) {
        return _mm_cmp_pd(a, b, i);
    }
    template<int i>
    static inline m256d vcmppd(m256d a, m256d b) {
        return _mm256_cmp_pd(a, b, i);
    }
    template<int i>
    static inline m128 vcmpss(m128 a, m128 b) {
        return _mm_cmp_ss(a, b, i);
    }
    template<int i>
    static inline m128d vcmpsd(m128d a, m128d b) {
        return _mm_cmp_sd(a, b, i);
    }

    static inline int vtestcps128(m128 a, m128 b) {
        return _mm_testc_ps(a, b);
    }
    static inline int vtestcps(m256 a, m256 b) {
        return _mm256_testc_ps(a, b);
    }
    static inline int vtestcpd128(m128d a, m128d b) {
        return _mm_testc_pd(a, b);
    }
    static inline int vtestcpd(m256d a, m256d b) {
        return _mm256_testc_pd(a, b);
    }
    static inline int vptestc(m256i a, m256i b) {
        return _mm256_testc_si256(a, b);
    }
    static inline int vtestnzcps128(m128 a, m128 b) {
        return _mm_testnzc_ps(a, b);
    }
    static inline int vtestnzcps(m256 a, m256 b) {
        return _mm256_testnzc_ps(a, b);
    }
    static inline int vtestnzcpd128(m128d a, m128d b) {
        return _mm_testnzc_pd(a, b);
    }
    static inline int vtestnzcpd(m256d a, m256d b) {
        return _mm256_testnzc_pd(a, b);
    }
    static inline int vptestnzc(m256i a, m256i b) {
        return _mm256_testnzc_si256(a, b);
    }
    static inline int vtestzps128(m128 a, m128 b) {
        return _mm_testz_ps(a, b);
    }
    static inline int vtestzps(m256 a, m256 b) {
        return _mm256_testz_ps(a, b);
    }
    static inline int vtestzpd128(m128d a, m128d b) {
        return _mm_testz_pd(a, b);
    }
    static inline int vtestzpd(m256d a, m256d b) {
        return _mm256_testz_pd(a, b);
    }
    static inline int vptestz(m256i a, m256i b) {
        return _mm256_testz_si256(a, b);
    }

    // maximum, minimum
    static inline m256 vmaxps(m256 a, m256 b){
        return _mm256_max_ps(a, b);
    }
    static inline m256d vmaxpd(m256d a, m256d b){
        return _mm256_max_pd(a, b);
    }
    static inline m256 vminps(m256 a, m256 b){
        return _mm256_min_ps(a, b);
    }
    static inline m256d vminpd(m256d a, m256d b){
        return _mm256_min_pd(a, b);
    }

}

#endif
