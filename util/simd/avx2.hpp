#ifndef SIMD_AVX2_HPP
#define SIMD_AVX2_HPP

#include "avx.hpp"

#include "immintrin.h"

namespace simd {

    // broadcast
    static inline m128i vpbroadcastb128(m128i a) {
        return _mm_broadcastb_epi8(a);
    }
    static inline m256i vpbroadcastb(m128i a) {
        return _mm256_broadcastb_epi8(a);
    }
    static inline m128i vpbroadcastw128(m128i a) {
        return _mm_broadcastw_epi16(a);
    }
    static inline m256i vpbroadcastw(m128i a) {
        return _mm256_broadcastw_epi16(a);
    }
    static inline m128i vpbroadcastd128(m128i a) {
        return _mm_broadcastd_epi32(a);
    }
    static inline m256i vpbroadcastd(m128i a) {
        return _mm256_broadcastd_epi32(a);
    }
    static inline m128i vpbroadcastq128(m128i a) {
        return _mm_broadcastq_epi64(a);
    }
    static inline m256i vpbroadcastq(m128i a) {
        return _mm256_broadcastq_epi64(a);
    }
    /*static inline m128d movddup(m128d a) {
        return _mm_broadcastsd_pd(a);
    }*/
    static inline m256d vbroadcastsd(m128d a) {
        return _mm256_broadcastsd_pd(a);
    }
    static inline m256i vbroadcasti128(m128i a) {
        return _mm256_broadcastsi128_si256(a);
    }
    static inline m128 vbroadcastss128(m128 a) {
        return _mm_broadcastss_ps(a);
    }
    static inline m256 vbroadcastss(m128 a) {
        return _mm256_broadcastss_ps(a);
    }

    // move and shuffle
    static inline m128i vpmaskmovd128(const int * a, m128i msk){
        return _mm_maskload_epi32(a, msk);
    }
    static inline m256i vpmaskmovd(const int * a, m256i msk) {
        return _mm256_maskload_epi32(a, msk);
    }
    static inline m128i vpmaskmovq128(const int64_t * a, m128i msk) {
        return _mm_maskload_epi64((const long long *)a, msk);
    }
    static inline m256i vpmaskmovq(const int64_t * a, m256i msk) {
        return _mm256_maskload_epi64((const long long *)a, msk);
    }
    static inline void vpmaskmovd128(int * a, m128i msk, m128i b) {
        _mm_maskstore_epi32(a, msk, b);
    }
    static inline void vpmaskmovd(int * a, m256i msk, m256i b) {
        _mm256_maskstore_epi32(a, msk, b);
    }
    static inline void vpmaskmovq128(int64_t * a, m128i msk, m128i b) {
        _mm_maskstore_epi64((long long *)a, msk, b);
    }
    static inline void vpmaskmovq(int64_t * a, m256i msk, m256i b) {
        _mm256_maskstore_epi64((long long *)a, msk, b);
    }
    static inline int vpmovmskb(m256i a) {
        return _mm256_movemask_epi8(a);
    }
    static inline m256i vmovntdqa(const m256i * a) {
        return _mm256_stream_load_si256(a);
    }
    template<int i>
    static inline m256i vperm2i128(m256i a, m256i b){
        return _mm256_permute2x128_si256(a, b, i);
    }
    template<int i>
    static inline m256i vpermq(m256i a) {
        return _mm256_permute4x64_epi64(a, i);
    }
    template<int i>
    static inline m256d vpermpd(m256d a) {
        return _mm256_permute4x64_pd(a, i);
    }
    static inline m256i vpermd(m256i a, m256i b) {
        return _mm256_permutevar8x32_epi32(a, b);
    }
    static inline m256 vpermps(m256 a, m256i b) {
        return _mm256_permutevar8x32_ps(a, b);
    }
    static inline m256i vpshufb(m256i a, m256i b) {
        return _mm256_shuffle_epi8(a, b);
    }
    template<int i>
    static inline m256i vpshufhw(m256i a) {
        return _mm256_shufflehi_epi16(a, i);
    }
    template<int i0, int i1, int i2, int i3>
    static inline m256i vpshufhw(m256i a) {
        return _mm256_shufflehi_epi16(a, i0 | (i1 << 2) | (i2 << 4) | (i3 << 6));
    }
    template<int i>
    static inline m256i vpshuflw(m256i a) {
        return _mm256_shufflelo_epi16(a, i);
    }
    template<int i0, int i1, int i2, int i3>
    static inline m256i vpshuflw(m256i a) {
        return _mm256_shufflelo_epi16(a, i0 | (i1 << 2) | (i2 << 4) | (i3 << 6));
    }
    template<int i>
    static inline m256i vpshufd(m256i a) {
        return _mm256_shuffle_epi32(a, i);
    }
    template<int i0, int i1, int i2, int i3>
    static inline m256i vpshufd(m256i a) {
        return _mm256_shuffle_epi32(a, i0 | (i1 << 2) | (i2 << 4) | (i3 << 6));
    }

    // pack
    static inline m256i vpacksswb(m256i a, m256i b) {
        return _mm256_packs_epi16(a, b);
    }
    static inline m256i vpackssdw(m256i a, m256i b) {
        return _mm256_packs_epi32(a, b);
    }
    static inline m256i vpackuswb(m256i a, m256i b) {
        return _mm256_packus_epi16(a, b);
    }
    static inline m256i vpackusd(m256i a, m256i b) {
        return _mm256_packus_epi32(a, b);
    }

    // unpack
    static inline m256i vpunpckhwd(m256i a, m256i b) {
        return _mm256_unpackhi_epi16(a, b);
    }
    static inline m256i vpunpckhdq(m256i a, m256i b) {
        return _mm256_unpackhi_epi32(a, b);
    }
    static inline m256i vpunpckhqdq(m256i a, m256i b) {
        return _mm256_unpackhi_epi64(a, b);
    }
    static inline m256i vpunpckhbw(m256i a, m256i b) {
        return _mm256_unpackhi_epi8(a, b);
    }
    static inline m256i vpunpcklwd(m256i a, m256i b) {
        return _mm256_unpacklo_epi16(a, b);
    }
    static inline m256i vpunpckldq(m256i a, m256i b) {
        return _mm256_unpacklo_epi32(a, b);
    }
    static inline m256i vpunpcklqdq(m256i a, m256i b) {
        return _mm256_unpacklo_epi64(a, b);
    }
    static inline m256i vpunpcklbw(m256i a, m256i b) {
        return _mm256_unpacklo_epi8(a, b);
    }

    // gather
    template<int i>
    static inline m128i vpgatherdd(const int * a, m128i b) {
        return _mm_i32gather_epi32(a, b, i);
    }
    template<int i>
    static inline m128i vpgatherdd(m128i a, const int * b, m128i c, m128i msk) {
        return _mm_mask_i32gather_epi32(a, b, c, msk, i);
    }
    template<int i>
    static inline m256i vpgatherdd(const int * a, m256i b) {
        return _mm256_i32gather_epi32(a, b, i);
    }
    template<int i>
    static inline m256i vpgatherdd(m256i a, const int * b, m256i c, m256i msk) {
        return _mm256_mask_i32gather_epi32(a, b, c, msk, i);
    }
    template<int i>
    static inline m128i vpgatherdq(const int64_t * a, m128i b) {
        return _mm_i32gather_epi64((long long *)a, b, i);
    }
    template<int i>
    static inline m128i vpgatherdq(m128i a, const int64_t * b, m128i c, m128i msk) {
        return _mm_mask_i32gather_epi64(a, (long long *)b, c, msk, i);
    }
    template<int i>
    static inline m256i vpgatherdq(const int64_t * a, m128i b) {
        return _mm256_i32gather_epi64((long long *)a, b, i);
    }
    template<int i>
    static inline m256i vpgatherdq(m256i a, const int64_t * b, m128i c, m256i msk) {
        return _mm256_mask_i32gather_epi64(a, b, c, msk, i);
    }
    template<int i>
    static inline m128d vgatherdpd(const double * a, m128i b) {
        return _mm_i32gather_pd(a, b, i);
    }
    template<int i>
    static inline m128d vgatherdpd(m128d a, const double * b, m128i c, m128d msk) {
        return _mm_mask_i32gather_pd(a, b, c, msk, i);
    }
    template<int i>
    static inline m256d vgatherdpd(const double * a, m128i b) {
        return _mm256_i32gather_pd(a, b, i);
    }
    template<int i>
    static inline m256d vgatherdp(m256d a, const double * b, m128i c, m256d msk) {
        return _mm256_mask_i32gather_pd(a, b, c, msk, i);
    }
    template<int i>
    static inline m128 vgatherdps(const float * a, m128i b) {
        return _mm_i32gather_ps(a, b, i);
    }
    template<int i>
    static inline m128 vgatherdps(m128 a, const float * b, m128i c, m128 msk) {
        return _mm_mask_i32gather_ps(a, b, c, msk, i);
    }
    template<int i>
    static inline m256 vgatherdps(const float * a, m256i b) {
        return _mm256_i32gather_ps(a, b, i);
    }
    template<int i>
    static inline m256 vgatherdps(m256 a, const float * b, m256i c, m256 msk) {
        return _mm256_mask_i32gather_ps(a, b, c, msk, i);
    }
    template<int i>
    static inline m128i vpgatherqd(const int * a, m128i b) {
        return _mm_i64gather_epi32(a, b, i);
    }
    template<int i>
    static inline m128i vpgatherqd(m128i a, const int * b, m128i c, m128i msk) {
        return _mm_mask_i64gather_epi32(a, b, c, msk, i);
    }
    template<int i>
    static inline m128i vpgatherqd(const int * a, m256i b) {
        return _mm256_i64gather_epi32(a, b, i);
    }
    template<int i>
    static inline m128i vpgatherqd(m128i a, const int * b, m256i c, m128i msk) {
        return _mm256_mask_i64gather_epi32(a, b, c, msk, i);
    }
    template<int i>
    static inline m128i vpgatherqq(const int64_t * a, m128i b) {
        return _mm_i64gather_epi64((long long *)a, b, i);
    }
    template<int i>
    static inline m128i vpgatherqq(m128i a, const int64_t * b, m128i c, m128i msk) {
        return _mm_mask_i64gather_epi64(a, (long long *)b, c, msk, i);
    }
    template<int i>
    static inline m256i vpgatherqq(const int64_t * a, m256i b) {
        return _mm256_i64gather_epi64((long long * )a, b, i);
    }
    template<int i>
    static inline m256i vpgatherqq(m256i a, const int64_t * b, m256i c, m256i msk) {
        return _mm256_mask_i64gather_epi64(a, (long long *)b, c, msk, i);
    }
    template<int i>
    static inline m128d vgatherqpd(const double * a, m128i b) {
        return _mm_i64gather_pd(a, b, i);
    }
    template<int i>
    static inline m128d vgatherqpd(m128d a, const double * b, m128i c, m128d msk) {
        return _mm_mask_i64gather_pd(a, b, c, msk, i);
    }
    template<int i>
    static inline m256d vgatherqpd(const double * a, m256i b) {
        return _mm256_i64gather_pd(a, b, i);
    }
    template<int i>
    static inline m256d vgatherqpd(m256d a, const double * b, m256i c, m256d msk) {
        return _mm256_mask_i64gather_pd(a, b, c, msk, i);
    }
    template<int i>
    static inline m128 vgatherqps(const float * a, m128i b) {
        return _mm_i64gather_ps(a, b, i);
    }
    template<int i>
    static inline m128 vgatherqps(m128 a, const float * b, m128i c, m128 msk) {
        return _mm_mask_i64gather_ps(a, b, c, msk, i);
    }
    template<int i>
    static inline m128 vgatherqps(const float * a, m256i b) {
        return _mm256_i64gather_ps(a, b, i);
    }
    template<int i>
    static inline m128 vgatherqps(m128 a, const float * b, m256i c, m128 msk) {
        return _mm256_mask_i64gather_ps(a, b, c, msk, i);
    }

    // align, blend, extract, insert
    template<int i>
    static inline m256i vpalignr(m256i a, m256i b) {
        return _mm256_alignr_epi8(a, b, i);
    }
    template<int i>
    static inline m256i vpblendw(m256i a, m256i b) {
        return _mm256_blend_epi16(a, b, i);
    }
    template<int i>
    static inline m128i vpblendd(m128i a, m128i b) {
        return _mm_blend_epi32(a, b, i);
    }
    template<int i>
    static inline m256i vpblendd(m256i a, m256i b) {
        return _mm256_blend_epi32(a, b, i);
    }
    static inline m256i vpblendvb(m256i a, m256i b, m256i msk) {
        return _mm256_blendv_epi8(a, b, msk);
    }
    template<int i>
    static inline m128i vextracti128(m256i a) {
        return _mm256_extracti128_si256(a, i);
    }
    template<int i>
    static inline m256i vinserti128(m256i a, m128i b) {
        return _mm256_inserti128_si256(a, b, i);
    }

    // addition
    static inline m256i vpaddb(m256i a, m256i b) {
        return _mm256_add_epi8(a, b);
    }
    static inline m256i vpaddw(m256i a, m256i b) {
        return _mm256_add_epi16(a, b);
    }
    static inline m256i vpaddd(m256i a, m256i b) {
        return _mm256_add_epi32(a, b);
    }
    static inline m256i vpaddq(m256i a, m256i b) {
        return _mm256_add_epi64(a, b);
    }
    static inline m256i vpaddsb(m256i a, m256i b) {
        return _mm256_adds_epi8(a, b);
    }
    static inline m256i vpaddsw(m256i a, m256i b) {
        return _mm256_adds_epi16(a, b);
    }
    static inline m256i vpaddusb(m256i a, m256i b) {
        return _mm256_adds_epu8(a, b);
    }
    static inline m256i vpaddusw(m256i a, m256i b) {
        return _mm256_adds_epu16(a, b);
    }
    static inline m256i vphaddw(m256i a, m256i b) {
        return _mm256_hadd_epi16(a, b);
    }
    static inline m256i vphaddd(m256i a, m256i b) {
        return _mm256_hadd_epi32(a, b);
    }
    static inline m256i vphaddsw(m256i a, m256i b) {
        return _mm256_hadds_epi16(a, b);
    }

    // subtraction
    static inline m256i vpsubw(m256i a, m256i b) {
        return _mm256_sub_epi16(a, b);
    }
    static inline m256i vpsubd(m256i a, m256i b) {
        return _mm256_sub_epi32(a, b);
    }
    static inline m256i vpsubq(m256i a, m256i b) {
        return _mm256_sub_epi64(a, b);
    }
    static inline m256i vpsubb(m256i a, m256i b) {
        return _mm256_sub_epi8(a, b);
    }
    static inline m256i vpsubsw(m256i a, m256i b) {
        return _mm256_subs_epi16(a, b);
    }
    static inline m256i vpsubsb(m256i a, m256i b) {
        return _mm256_subs_epi8(a, b);
    }
    static inline m256i vpsubusw(m256i a, m256i b) {
        return _mm256_subs_epu16(a, b);
    }
    static inline m256i vpsubusb(m256i a, m256i b) {
        return _mm256_subs_epu8(a, b);
    }
    static inline m256i vphsubw(m256i a, m256i b) {
        return _mm256_hsub_epi16(a, b);
    }
    static inline m256i vphsubd(m256i a, m256i b) {
        return _mm256_hsub_epi32(a, b);
    }
    static inline m256i vphsubsw(m256i a, m256i b) {
        return _mm256_hsubs_epi16(a, b);
    }
    template<int i>
    static inline m256i vmpsadbw(m256i a, m256i b) {
        return _mm256_mpsadbw_epu8(a, b, i);
    }
    static inline m256i vpsadbw(m256i a, m256i b) {
        return _mm256_sad_epu8(a, b);
    }
    static inline m256i vpsignw(m256i a, m256i b) {
        return _mm256_sign_epi16(a, b);
    }
    static inline m256i vpsignd(m256i a, m256i b) {
        return _mm256_sign_epi32(a, b);
    }
    static inline m256i vpsignb(m256i a, m256i b) {
        return _mm256_sign_epi8(a, b);
    }

    // multiplication
    static inline m256i vpmuldq(m256i a, m256i b) {
        return _mm256_mul_epi32(a, b);
    }
    static inline m256i vpmuludq(m256i a, m256i b) {
        return _mm256_mul_epu32(a, b);
    }
    static inline m256i vpmulhw(m256i a, m256i b) {
        return _mm256_mulhi_epi16(a, b);
    }
    static inline m256i vpmulhuw(m256i a, m256i b) {
        return _mm256_mulhi_epu16(a, b);
    }
    static inline m256i vpmulhrsw(m256i a, m256i b) {
        return _mm256_mulhrs_epi16(a, b);
    }
    static inline m256i vpmullw(m256i a, m256i b) {
        return _mm256_mullo_epi16(a, b);
    }
    static inline m256i vpmulld(m256i a, m256i b) {
        return _mm256_mullo_epi32(a, b);
    }
    static inline m256i vpmaddwd(m256i a, m256i b) {
        return _mm256_madd_epi16(a, b);
    }
    static inline m256i vpmaddubsw(m256i a, m256i b) {
        return _mm256_maddubs_epi16(a, b);
    }

    // average
    static inline m256i vpavgw(m256i a, m256i b) {
        return _mm256_avg_epu16(a, b);
    }
    static inline m256i vpavgb(m256i a, m256i b) {
        return _mm256_avg_epu8(a, b);
    }

    // absolute value
    static inline m256i vpabsw(m256i a) {
        return _mm256_abs_epi16(a);
    }
    static inline m256i vpabsd(m256i a) {
        return _mm256_abs_epi32(a);
    }
    static inline m256i vpabsb(m256i a) {
        return _mm256_abs_epi8(a);
    }

    // logic
    static inline m256i vpand(m256i a, m256i b) {
        return _mm256_and_si256(a, b);
    }
    static inline m256i vpandn(m256i a, m256i b) {
        return _mm256_andnot_si256(a, b);
    }
    static inline m256i vpor(m256i a, m256i b) {
        return _mm256_or_si256(a, b);
    }
    static inline m256i vpxor(m256i a, m256i b) {
        return _mm256_xor_si256(a, b);
    }

    // bit shift
    template<int i>
    static inline m256i vpslldq(m256i a) {
        return _mm256_bslli_epi128(a, i);
    }
    template<int i>
    static inline m256i vpsrldq(m256i a) {
        return _mm256_bsrli_epi128(a, i);
    }
    static inline m256i vpsllw(m256i a, m128i b) {
        return _mm256_sll_epi16(a, b);
    }
    static inline m256i vpslld(m256i a, m128i b) {
        return _mm256_sll_epi32(a, b);
    }
    static inline m256i vpsllq(m256i a, m128i b) {
        return _mm256_sll_epi64(a, b);
    }
    template<int i>
    static inline m256i vpsllw(m256i a) {
        return _mm256_slli_epi16(a, i);
    }
    template<int i>
    static inline m256i vpslld(m256i a) {
        return _mm256_slli_epi32(a, i);
    }
    template<int i>
    static inline m256i vpsllq(m256i a) {
        return _mm256_slli_epi64(a, i);
    }
    static inline m128i vpsllvd128(m128i a, m128i b) {
        return _mm_sllv_epi32(a, b);
    }
    static inline m256i vpsllvd(m256i a, m256i b) {
        return _mm256_sllv_epi32(a, b);
    }
    static inline m128i vpsllvq128(m128i a, m128i b) {
        return _mm_sllv_epi64(a, b);
    }
    static inline m256i vpsllvq(m256i a, m256i b) {
        return _mm256_sllv_epi64(a, b);
    }
    static inline m256i vpsraw(m256i a, m128i b) {
        return _mm256_sra_epi16(a, b);
    }
    static inline m256i vpsrad(m256i a, m128i b) {
        return _mm256_sra_epi32(a, b);
    }
    template<int i>
    static inline m256i vpsraw(m256i a) {
        return _mm256_srai_epi16(a, i);
    }
    template<int i>
    static inline m256i vpsrad(m256i a) {
        return _mm256_srai_epi32(a, i);
    }
    static inline m128i vpsravd128(m128i a, m128i b) {
        return _mm_srav_epi32(a, b);
    }
    static inline m256i vpsravd(m256i a, m256i b) {
        return _mm256_srav_epi32(a, b);
    }
    static inline m256i vpsrlw(m256i a, m128i b) {
        return _mm256_srl_epi16(a, b);
    }
    static inline m256i vpsrld(m256i a, m128i b) {
        return _mm256_srl_epi32(a, b);
    }
    static inline m256i vpsrlq(m256i a, m128i b) {
        return _mm256_srl_epi64(a, b);
    }
    template<int i>
    static inline m256i vpsrlw(m256i a) {
        return _mm256_srli_epi16(a, i);
    }
    template<int i>
    static inline m256i vpsrld(m256i a) {
        return _mm256_srli_epi32(a, i);
    }
    template<int i>
    static inline m256i vpsrlq(m256i a) {
        return _mm256_srli_epi64(a, i);
    }
    static inline m128i vpsrlvd128(m128i a, m128i b) {
        return _mm_srlv_epi32(a, b);
    }
    static inline m256i vpsrlvd(m256i a, m256i b) {
        return _mm256_srlv_epi32(a, b);
    }
    static inline m128i vpsrlvq128(m128i a, m128i b) {
        return _mm_srlv_epi64(a, b);
    }
    static inline m256i vpsrlvq(m256i a, m256i b) {
        return _mm256_srlv_epi64(a, b);
    }

    // sign and zero extend
    static inline m256i vpmovsxwd(m128i a) {
        return _mm256_cvtepi16_epi32(a);
    }
    static inline m256i vpmovsxwq(m128i a) {
        return _mm256_cvtepi16_epi64(a);
    }
    static inline m256i vpmovsxdq(m128i a) {
        return _mm256_cvtepi32_epi64(a);
    }
    static inline m256i vpmovsxbw(m128i a) {
        return _mm256_cvtepi8_epi16(a);
    }
    static inline m256i vpmovsxbd(m128i a) {
        return _mm256_cvtepi8_epi32(a);
    }
    static inline m256i vpmovsxbq(m128i a) {
        return _mm256_cvtepi8_epi64(a);
    }
    static inline m256i vpmovzxwd(m128i a) {
        return _mm256_cvtepu16_epi32(a);
    }
    static inline m256i vpmovzxwq(m128i a) {
        return _mm256_cvtepu16_epi64(a);
    }
    static inline m256i vpmovzxdq(m128i a) {
        return _mm256_cvtepu32_epi64(a);
    }
    static inline m256i vpmovzxbw(m128i a) {
        return _mm256_cvtepu8_epi16(a);
    }
    static inline m256i vpmovzxbd(m128i a) {
        return _mm256_cvtepu8_epi32(a);
    }
    static inline m256i vpmovzxbq(m128i a) {
        return _mm256_cvtepu8_epi64(a);
    }

    // comparison
    static inline m256i vpcmpeqw(m256i a, m256i b) {
        return _mm256_cmpeq_epi16(a, b);
    }
    static inline m256i vpcmpeqd(m256i a, m256i b) {
        return _mm256_cmpeq_epi32(a, b);
    }
    static inline m256i vpcmpeqq(m256i a, m256i b) {
        return _mm256_cmpeq_epi64(a, b);
    }
    static inline m256i vpcmpeqb(m256i a, m256i b) {
        return _mm256_cmpeq_epi8(a, b);
    }
    static inline m256i vpcmpgtw(m256i a, m256i b) {
        return _mm256_cmpgt_epi16(a, b);
    }
    static inline m256i vpcmpgtd(m256i a, m256i b) {
        return _mm256_cmpgt_epi32(a, b);
    }
    static inline m256i vpcmpgtq(m256i a, m256i b) {
        return _mm256_cmpgt_epi64(a, b);
    }
    static inline m256i vpcmpgtb(m256i a, m256i b) {
        return _mm256_cmpgt_epi8(a, b);
    }

    // maximum, minimum
    static inline m256i vpmaxsw(m256i a, m256i b) {
        return _mm256_max_epi16(a, b);
    }
    static inline m256i vpmaxsd(m256i a, m256i b) {
        return _mm256_max_epi32(a, b);
    }
    static inline m256i vpmaxsb(m256i a, m256i b) {
        return _mm256_max_epi8(a, b);
    }
    static inline m256i vpmaxuw(m256i a, m256i b) {
        return _mm256_max_epu16(a, b);
    }
    static inline m256i vpmaxud(m256i a, m256i b) {
        return _mm256_max_epu32(a, b);
    }
    static inline m256i vpmaxub(m256i a, m256i b) {
        return _mm256_max_epu8(a, b);
    }
    static inline m256i vpminsw(m256i a, m256i b) {
        return _mm256_min_epi16(a, b);
    }
    static inline m256i vpminsd(m256i a, m256i b) {
        return _mm256_min_epi32(a, b);
    }
    static inline m256i vpminsb(m256i a, m256i b) {
        return _mm256_min_epi8(a, b);
    }
    static inline m256i vpminuw(m256i a, m256i b) {
        return _mm256_min_epu16(a, b);
    }
    static inline m256i vpminud(m256i a, m256i b) {
        return _mm256_min_epu32(a, b);
    }
    static inline m256i vpminub(m256i a, m256i b) {
        return _mm256_min_epu8(a, b);
    }

}

#endif
