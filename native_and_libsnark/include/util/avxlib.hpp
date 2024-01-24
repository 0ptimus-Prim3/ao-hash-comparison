#pragma once

#include <cinttypes>
#include <iostream>

#include "intrinsics.h"

inline std::ostream &operator<<(std::ostream &os, __m128i v)
{
    uint8_t a[sizeof(v)]{};

    _mm_storeu_si128((__m128i *)a, v);

    os << '{';
    for (size_t i = 0; i < sizeof(v) - 1; ++i)
        os << (int)a[(sizeof(v) - 1) - i] << ", ";

    return os << (int)a[0] << '}';
}

inline std::ostream &operator<<(std::ostream &os, __m256i v)
{
    uint8_t a[sizeof(v)]{};

    _mm256_storeu_si256((__m256i *)a, v);

    os << '{';
    for (size_t i = 0; i < sizeof(v) - 1; ++i)
        os << (int)a[(sizeof(v) - 1) - i] << ", ";

    return os << (int)a[0] << '}';
}

inline std::ostream &operator<<(std::ostream &os, __m512i v)
{
    uint8_t a[sizeof(v)]{};

    _mm512_storeu_si512((__m512i *)a, v);

    os << '{';
    for (size_t i = 0; i < sizeof(v) - 1; ++i)
        os << (int)a[(sizeof(v) - 1) - i] << ", ";

    return os << (int)a[0] << '}';
}

namespace avx
{
    // clang-format off
#define mm256_ptable_to_pbytes(t) _mm256_set_epi8(                                      \
    (char)(t[31] == 0xFF ? 0xFF : t[31] / 8), (char)(t[30] == 0xFF ? 0xFF : t[30] / 8), \
    (char)(t[29] == 0xFF ? 0xFF : t[29] / 8), (char)(t[28] == 0xFF ? 0xFF : t[28] / 8), \
    (char)(t[27] == 0xFF ? 0xFF : t[27] / 8), (char)(t[26] == 0xFF ? 0xFF : t[26] / 8), \
    (char)(t[25] == 0xFF ? 0xFF : t[25] / 8), (char)(t[24] == 0xFF ? 0xFF : t[24] / 8), \
    (char)(t[23] == 0xFF ? 0xFF : t[23] / 8), (char)(t[22] == 0xFF ? 0xFF : t[22] / 8), \
    (char)(t[21] == 0xFF ? 0xFF : t[21] / 8), (char)(t[20] == 0xFF ? 0xFF : t[20] / 8), \
    (char)(t[19] == 0xFF ? 0xFF : t[19] / 8), (char)(t[18] == 0xFF ? 0xFF : t[18] / 8), \
    (char)(t[17] == 0xFF ? 0xFF : t[17] / 8), (char)(t[16] == 0xFF ? 0xFF : t[16] / 8), \
    (char)(t[15] == 0xFF ? 0xFF : t[15] / 8), (char)(t[14] == 0xFF ? 0xFF : t[14] / 8), \
    (char)(t[13] == 0xFF ? 0xFF : t[13] / 8), (char)(t[12] == 0xFF ? 0xFF : t[12] / 8), \
    (char)(t[11] == 0xFF ? 0xFF : t[11] / 8), (char)(t[10] == 0xFF ? 0xFF : t[10] / 8), \
    (char)(t[9]  == 0xFF ? 0xFF : t[9]  / 8), (char)(t[8]  == 0xFF ? 0xFF : t[8]  / 8), \
    (char)(t[7]  == 0xFF ? 0xFF : t[7]  / 8), (char)(t[6]  == 0xFF ? 0xFF : t[6]  / 8), \
    (char)(t[5]  == 0xFF ? 0xFF : t[5]  / 8), (char)(t[4]  == 0xFF ? 0xFF : t[4]  / 8), \
    (char)(t[3]  == 0xFF ? 0xFF : t[3]  / 8), (char)(t[2]  == 0xFF ? 0xFF : t[2]  / 8), \
    (char)(t[1]  == 0xFF ? 0xFF : t[1]  / 8), (char)(t[0]  == 0xFF ? 0xFF : t[0]  / 8))

#define mm256_ptable_to_pbits(t) _mm256_set_epi8(       \
    (char)(1 << (t[31] % 8)), (char)(1 << (t[30] % 8)), \
    (char)(1 << (t[29] % 8)), (char)(1 << (t[28] % 8)), \
    (char)(1 << (t[27] % 8)), (char)(1 << (t[26] % 8)), \
    (char)(1 << (t[25] % 8)), (char)(1 << (t[24] % 8)), \
    (char)(1 << (t[23] % 8)), (char)(1 << (t[22] % 8)), \
    (char)(1 << (t[21] % 8)), (char)(1 << (t[20] % 8)), \
    (char)(1 << (t[19] % 8)), (char)(1 << (t[18] % 8)), \
    (char)(1 << (t[17] % 8)), (char)(1 << (t[16] % 8)), \
    (char)(1 << (t[15] % 8)), (char)(1 << (t[14] % 8)), \
    (char)(1 << (t[13] % 8)), (char)(1 << (t[12] % 8)), \
    (char)(1 << (t[11] % 8)), (char)(1 << (t[10] % 8)), \
    (char)(1 << (t[9]  % 8)), (char)(1 << (t[8]  % 8)), \
    (char)(1 << (t[7]  % 8)), (char)(1 << (t[6]  % 8)), \
    (char)(1 << (t[5]  % 8)), (char)(1 << (t[4]  % 8)), \
    (char)(1 << (t[3]  % 8)), (char)(1 << (t[2]  % 8)), \
    (char)(1 << (t[1]  % 8)), (char)(1 << (t[0]  % 8)))

#define mm512_ptable_to_pbytes(t) _mm512_set_epi8(                                      \
    (char)(t[63] == 0xFF ? 0xFF : t[63] / 8), (char)(t[62] == 0xFF ? 0xFF : t[62] / 8), \
    (char)(t[61] == 0xFF ? 0xFF : t[61] / 8), (char)(t[60] == 0xFF ? 0xFF : t[60] / 8), \
    (char)(t[59] == 0xFF ? 0xFF : t[59] / 8), (char)(t[58] == 0xFF ? 0xFF : t[58] / 8), \
    (char)(t[57] == 0xFF ? 0xFF : t[57] / 8), (char)(t[56] == 0xFF ? 0xFF : t[56] / 8), \
    (char)(t[55] == 0xFF ? 0xFF : t[55] / 8), (char)(t[54] == 0xFF ? 0xFF : t[54] / 8), \
    (char)(t[53] == 0xFF ? 0xFF : t[53] / 8), (char)(t[52] == 0xFF ? 0xFF : t[52] / 8), \
    (char)(t[51] == 0xFF ? 0xFF : t[51] / 8), (char)(t[50] == 0xFF ? 0xFF : t[50] / 8), \
    (char)(t[49] == 0xFF ? 0xFF : t[49] / 8), (char)(t[48] == 0xFF ? 0xFF : t[48] / 8), \
    (char)(t[47] == 0xFF ? 0xFF : t[47] / 8), (char)(t[46] == 0xFF ? 0xFF : t[46] / 8), \
    (char)(t[45] == 0xFF ? 0xFF : t[45] / 8), (char)(t[44] == 0xFF ? 0xFF : t[44] / 8), \
    (char)(t[43] == 0xFF ? 0xFF : t[43] / 8), (char)(t[42] == 0xFF ? 0xFF : t[42] / 8), \
    (char)(t[41] == 0xFF ? 0xFF : t[41] / 8), (char)(t[40] == 0xFF ? 0xFF : t[40] / 8), \
    (char)(t[39] == 0xFF ? 0xFF : t[39] / 8), (char)(t[38] == 0xFF ? 0xFF : t[38] / 8), \
    (char)(t[37] == 0xFF ? 0xFF : t[37] / 8), (char)(t[36] == 0xFF ? 0xFF : t[36] / 8), \
    (char)(t[35] == 0xFF ? 0xFF : t[35] / 8), (char)(t[34] == 0xFF ? 0xFF : t[34] / 8), \
    (char)(t[33] == 0xFF ? 0xFF : t[33] / 8), (char)(t[32] == 0xFF ? 0xFF : t[32] / 8), \
    (char)(t[31] == 0xFF ? 0xFF : t[31] / 8), (char)(t[30] == 0xFF ? 0xFF : t[30] / 8), \
    (char)(t[29] == 0xFF ? 0xFF : t[29] / 8), (char)(t[28] == 0xFF ? 0xFF : t[28] / 8), \
    (char)(t[27] == 0xFF ? 0xFF : t[27] / 8), (char)(t[26] == 0xFF ? 0xFF : t[26] / 8), \
    (char)(t[25] == 0xFF ? 0xFF : t[25] / 8), (char)(t[24] == 0xFF ? 0xFF : t[24] / 8), \
    (char)(t[23] == 0xFF ? 0xFF : t[23] / 8), (char)(t[22] == 0xFF ? 0xFF : t[22] / 8), \
    (char)(t[21] == 0xFF ? 0xFF : t[21] / 8), (char)(t[20] == 0xFF ? 0xFF : t[20] / 8), \
    (char)(t[19] == 0xFF ? 0xFF : t[19] / 8), (char)(t[18] == 0xFF ? 0xFF : t[18] / 8), \
    (char)(t[17] == 0xFF ? 0xFF : t[17] / 8), (char)(t[16] == 0xFF ? 0xFF : t[16] / 8), \
    (char)(t[15] == 0xFF ? 0xFF : t[15] / 8), (char)(t[14] == 0xFF ? 0xFF : t[14] / 8), \
    (char)(t[13] == 0xFF ? 0xFF : t[13] / 8), (char)(t[12] == 0xFF ? 0xFF : t[12] / 8), \
    (char)(t[11] == 0xFF ? 0xFF : t[11] / 8), (char)(t[10] == 0xFF ? 0xFF : t[10] / 8), \
    (char)(t[9]  == 0xFF ? 0xFF : t[9]  / 8), (char)(t[8]  == 0xFF ? 0xFF : t[8]  / 8), \
    (char)(t[7]  == 0xFF ? 0xFF : t[7]  / 8), (char)(t[6]  == 0xFF ? 0xFF : t[6]  / 8), \
    (char)(t[5]  == 0xFF ? 0xFF : t[5]  / 8), (char)(t[4]  == 0xFF ? 0xFF : t[4]  / 8), \
    (char)(t[3]  == 0xFF ? 0xFF : t[3]  / 8), (char)(t[2]  == 0xFF ? 0xFF : t[2]  / 8), \
    (char)(t[1]  == 0xFF ? 0xFF : t[1]  / 8), (char)(t[0]  == 0xFF ? 0xFF : t[0]  / 8))

#define mm512_ptable_to_pbits(t) _mm512_set_epi8(       \
    (char)(1 << (t[63] % 8)), (char)(1 << (t[62] % 8)), \
    (char)(1 << (t[61] % 8)), (char)(1 << (t[60] % 8)), \
    (char)(1 << (t[59] % 8)), (char)(1 << (t[58] % 8)), \
    (char)(1 << (t[57] % 8)), (char)(1 << (t[56] % 8)), \
    (char)(1 << (t[55] % 8)), (char)(1 << (t[54] % 8)), \
    (char)(1 << (t[53] % 8)), (char)(1 << (t[52] % 8)), \
    (char)(1 << (t[51] % 8)), (char)(1 << (t[50] % 8)), \
    (char)(1 << (t[49] % 8)), (char)(1 << (t[48] % 8)), \
    (char)(1 << (t[47] % 8)), (char)(1 << (t[46] % 8)), \
    (char)(1 << (t[45] % 8)), (char)(1 << (t[44] % 8)), \
    (char)(1 << (t[43] % 8)), (char)(1 << (t[42] % 8)), \
    (char)(1 << (t[41] % 8)), (char)(1 << (t[40] % 8)), \
    (char)(1 << (t[39] % 8)), (char)(1 << (t[38] % 8)), \
    (char)(1 << (t[37] % 8)), (char)(1 << (t[36] % 8)), \
    (char)(1 << (t[35] % 8)), (char)(1 << (t[34] % 8)), \
    (char)(1 << (t[33] % 8)), (char)(1 << (t[32] % 8)), \
    (char)(1 << (t[31] % 8)), (char)(1 << (t[30] % 8)), \
    (char)(1 << (t[29] % 8)), (char)(1 << (t[28] % 8)), \
    (char)(1 << (t[27] % 8)), (char)(1 << (t[26] % 8)), \
    (char)(1 << (t[25] % 8)), (char)(1 << (t[24] % 8)), \
    (char)(1 << (t[23] % 8)), (char)(1 << (t[22] % 8)), \
    (char)(1 << (t[21] % 8)), (char)(1 << (t[20] % 8)), \
    (char)(1 << (t[19] % 8)), (char)(1 << (t[18] % 8)), \
    (char)(1 << (t[17] % 8)), (char)(1 << (t[16] % 8)), \
    (char)(1 << (t[15] % 8)), (char)(1 << (t[14] % 8)), \
    (char)(1 << (t[13] % 8)), (char)(1 << (t[12] % 8)), \
    (char)(1 << (t[11] % 8)), (char)(1 << (t[10] % 8)), \
    (char)(1 << (t[9]  % 8)), (char)(1 << (t[8]  % 8)), \
    (char)(1 << (t[7]  % 8)), (char)(1 << (t[6]  % 8)), \
    (char)(1 << (t[5]  % 8)), (char)(1 << (t[4]  % 8)), \
    (char)(1 << (t[3]  % 8)), (char)(1 << (t[2]  % 8)), \
    (char)(1 << (t[1]  % 8)), (char)(1 << (t[0]  % 8)))


    // clang-format on

    static const __m128 zero_128 = _mm_setzero_ps();
    static const __m128d zero_128d = _mm_setzero_pd();
    static const __m128i zero_128i = _mm_setzero_si128();

    static const __m256 zero_256 = _mm256_setzero_ps();
    static const __m256d zero_256d = _mm256_setzero_pd();
    static const __m256i zero_256i = _mm256_setzero_si256();

#ifdef __AVX512F__
    static const __m512 zero_512 = _mm512_setzero_ps();
    static const __m512d zero_512d = _mm512_setzero_pd();
    static const __m512i zero_512i = _mm512_setzero_si512();
#endif

    // convert the bits of a type to the bytes in a wide enough register
    inline __m128i btob_no_norm(uint16_t x); // the bytes are 0/FF instead of 0/1
    inline __m256i btob_no_norm(uint32_t x);
    inline __m512i btob_no_norm(uint64_t x);
    inline uint64_t btob(uint8_t x);
    inline __m128i btob(uint16_t x);
    inline __m256i btob(uint32_t x);
    inline __m512i btob(uint64_t x);

    // permute the bits of x according to p. 64bit would require AVX512
    inline uint16_t permute(uint16_t x, __m128i pbits);
    inline uint32_t permute(uint32_t x, __m256i pbits, __m256i pbytes);
    inline uint64_t permute(uint64_t x, __m512i pbits, __m512i pbytes);

    // substitute the nibbles of x according to s.
    inline uint16_t sub_nib(uint16_t x, __m128i s);
    inline uint32_t sub_nib(uint32_t x, __m128i s);
    inline uint64_t sub_nib(uint64_t x, __m128i s);

    inline __m128i btob_no_norm(uint16_t x)
    {
        __m128i v_x = _mm_set1_epi16(x);
        const __m128i v_shuf = _mm_set_epi8(1, 1, 1, 1, 1, 1, 1, 1, //
                                            0, 0, 0, 0, 0, 0, 0, 0);
        __m128i v_and = _mm_set_epi8(0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01, //
                                     0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01);

        v_x = _mm_shuffle_epi8(v_x, v_shuf); // duplicate x[i] in x_v[8i:8i+7]
        v_x = _mm_and_si128(v_and, v_x);     // select the right bit for every byte
        v_x = _mm_cmpeq_epi8(v_and, v_x);    // normalize to 0/FF

        return v_x;
    }

    inline __m256i btob_no_norm(uint32_t x)
    {
        __m256i v_x = _mm256_set1_epi32(x);
        const __m256i v_shuf = _mm256_set_epi64x(0x0303030303030303, 0x0202020202020202, //
                                                 0x0101010101010101, 0x0000000000000000);
        __m256i v_and = _mm256_set1_epi64x(0x8040201008040201);

        v_x = _mm256_shuffle_epi8(v_x, v_shuf); // duplicate x[i] in x_v[8i:8i+7]
        v_x = _mm256_and_si256(v_and, v_x);     // select the right bit for every byte
        v_x = _mm256_cmpeq_epi8(v_and, v_x);    // normalize to 0/FF

        return v_x;
    }

    inline uint64_t btob(uint8_t x)
    {
        return _pdep_u64(x, 0x0101010101010101);
    }

    inline __m128i btob(uint16_t x)
    {
        return _mm_and_si128(_mm_set1_epi8(1), btob_no_norm(x)); // normalize to 0/1
    }

    inline __m256i btob(uint32_t x)
    {
        return _mm256_and_si256(_mm256_set1_epi8(1), btob_no_norm(x)); // normalize to 0/1
    }

    inline uint16_t permute(uint16_t x, __m128i pbits)
    {
        uint8_t l = (uint8_t)(x >> 8);
        uint8_t r = (uint8_t)x;
        __m128i v_x = _mm_set_epi8(l, l, l, l, l, l, l, l, r, r, r, r, r, r, r, r);
        __m128i msk = _mm_set_epi8((unsigned char)0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1, //
                                   (unsigned char)0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1);

        v_x = _mm_and_si128(v_x, msk);                  // select the right bit for every byte
        v_x = _mm_shuffle_epi8(v_x, pbits);             // move the bytes in the correct positions
        v_x = _mm_cmpeq_epi8(v_x, _mm_setzero_si128()); // set highest bit
        x = (uint16_t)~_mm_movemask_epi8(v_x);          // gather the result

        return x;
    }

    inline uint32_t permute(uint32_t x, __m256i pbits, __m256i pbytes)
    {
        __m256i v_x = _mm256_set1_epi32(x);

        v_x = _mm256_shuffle_epi8(v_x, pbytes); // move the bytes in the correct positions
        v_x = _mm256_and_si256(v_x, pbits);     // select the right bit for every byte
        v_x = _mm256_cmpeq_epi8(v_x, pbits);    // set highest bit
        x = _mm256_movemask_epi8(v_x);          // gather result

        return x;
    }

    inline uint64_t permute(uint64_t x, __m512i pbits, __m512i pbytes)
    {
        __m512i v_x = _mm512_set1_epi64(x);

        v_x = _mm512_shuffle_epi8(v_x, pbytes); // move the bytes in the correct positions
        v_x = _mm512_and_si512(v_x, pbits);     // select the right bit for every byte
        x = _mm512_cmpeq_epi8_mask(v_x, pbits); // set highest bit

        return x;
    }

    inline uint16_t sub_nib(uint16_t x, __m128i s)
    {
        __m128i v_x = _mm_cvtsi32_si128(_pdep_u32(x, 0x0F0F0F0F));

        v_x = _mm_shuffle_epi8(s, v_x);
        x = (uint16_t)_pext_u32(_mm_cvtsi128_si32(v_x), 0x0F0F0F0F);

        return x;
    }

    inline uint32_t sub_nib(uint32_t x, __m128i s)
    {
        __m128i v_x = _mm_cvtsi64_si128(_pdep_u64(x, 0x0F0F0F0F0F0F0F0F));

        v_x = _mm_shuffle_epi8(s, v_x);
        x = (uint32_t)_pext_u64(_mm_cvtsi128_si64(v_x), 0x0F0F0F0F0F0F0F0F);

        return x;
    }

    inline uint64_t sub_nib(uint64_t x, __m128i s)
    {
        uint64_t l = _pdep_u64(x >> 32, 0x0F0F0F0F0F0F0F0F);
        uint64_t r = _pdep_u64(x, 0x0F0F0F0F0F0F0F0F);
        __m128i v_x = _mm_set_epi64x(l, r);

        v_x = _mm_shuffle_epi8(s, v_x);
        l = _mm_extract_epi64(v_x, 1);
        r = _mm_cvtsi128_si64(v_x);
        x = _pext_u64(l, 0x0F0F0F0F0F0F0F0F) << 32 | _pext_u64(r, 0x0F0F0F0F0F0F0F0F);

        return x;
    }
}; // namespace avx
