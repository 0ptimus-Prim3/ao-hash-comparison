#pragma once

#include <array>
#include <cstddef>
#include <gmpxx.h>
#include <libff/common/default_types/ec_pp.hpp>
#include <utility>

enum class FieldKind
{
    UNKNOWN,
    LIBFF,
};

template<typename FieldT,
         std::enable_if_t<std::is_same_v<std::decay_t<FieldT>,
                                         libff::Fp_model<std::decay_t<FieldT>::num_limbs,
                                                         std::decay_t<FieldT>::mod>>> * = nullptr>
constexpr FieldKind field_kind()
{
    return FieldKind::LIBFF;
}

template<typename FieldT>
constexpr const char *field_kind_name()
{
    constexpr FieldKind kind = field_kind<FieldT>();

    if constexpr (kind == FieldKind::LIBFF)
        return "libff";
    else
        static_assert(false, "Invalid field type");
}

template<typename FieldT>
constexpr void field_init()
{
    // We initialize everything, tying NTL's (runtime) state to libff's (compile-time) state.
    libff::default_ec_pp::init_public_params();
    []
    {
        std::stringstream ss;

        ss << libff::Fr<libff::default_ec_pp>::mod;
    }();
}

template<typename FieldT>
constexpr size_t field_size()
{
    constexpr FieldKind kind = field_kind<FieldT>();

    if constexpr (kind == FieldKind::LIBFF)
        return FieldT::num_limbs * sizeof(mp_limb_t);
}

template<typename FieldT>
FieldT field_load(const void *src)
{
    static constexpr size_t FIELD_SIZE = field_size<FieldT>();
    static constexpr FieldKind kind = field_kind<FieldT>();

    if constexpr (kind == FieldKind::LIBFF)
    {
        mpz_class tmp;

        mpz_import(tmp.get_mpz_t(), FIELD_SIZE / sizeof(mp_limb_t), -1, sizeof(mp_limb_t), 0, 0,
                   src);

        return FieldT{tmp.get_mpz_t()};
    }
    else
        static_assert(false, "Invalid field type");
}

template<typename FieldT>
void field_store(void *dst, const FieldT &src)
{
    static constexpr size_t FIELD_SIZE = field_size<FieldT>();
    static constexpr FieldKind kind = field_kind<FieldT>();

    memset(dst, 0, FIELD_SIZE);

    if constexpr (kind == FieldKind::LIBFF)
    {
        mpz_class tmp;

        src.as_bigint().to_mpz(tmp.get_mpz_t());
        mpz_export(dst, NULL, -1, sizeof(mp_limb_t), 0, 0, tmp.get_mpz_t());
    }
    else
        static_assert(false, "Invalid field type");
}

template<typename Field>
void field_load(Field *dst, const void *src, size_t n)
{
    static constexpr size_t FIELD_SIZE = field_size<Field>();
    const unsigned char *data = reinterpret_cast<const unsigned char *>(src);

    for (size_t i = 0; i < n; ++i)
        dst[i] = field_load<Field>(data + FIELD_SIZE * i);
}

template<typename FieldT>
void field_store(void *dst, const FieldT *src, size_t n)
{
    static constexpr size_t FIELD_SIZE = field_size<FieldT>();
    unsigned char *data = reinterpret_cast<unsigned char *>(dst);

    memset(data, 0, FIELD_SIZE * n);

    for (size_t i = 0; i < n; ++i)
        field_store(data + FIELD_SIZE * i, src[i]);
}

template<typename Field>
Field field_inverse(const Field &x)
{
    constexpr FieldKind kind = field_kind<Field>();

    if constexpr (kind == FieldKind::LIBFF)
        return x.inverse();
}

template<typename Field>
constexpr Field field_random()
{
    constexpr FieldKind kind = field_kind<Field>();

    if constexpr (kind == FieldKind::LIBFF)
        return Field::random_element();
}

template<typename FieldT>
void field_clamp(void *vdata, size_t sz)
{
    static constexpr size_t FIELD_SIZE = field_size<FieldT>();

    mpz_class tmp;
    uint8_t *data = (uint8_t *)(vdata);

    for (size_t i = 0; i < sz / FIELD_SIZE; ++i)
    {
        mpz_import(tmp.get_mpz_t(), FIELD_SIZE, 1, 1, 0, 0, data + i * FIELD_SIZE);

        FieldT t{tmp.get_mpz_t()};
        t.as_bigint().to_mpz(tmp.get_mpz_t());

        memset(data + i * FIELD_SIZE, 0, FIELD_SIZE);
        mpz_export(data + i * FIELD_SIZE, NULL, 1, 1, 0, 0, tmp.get_mpz_t());
    }
}

template<typename FieldT>
mpz_class field_to_mpz(const FieldT &x)
{
    mpz_class x_mpz;

    x.as_bigint().to_mpz(x_mpz.get_mpz_t());

    return x_mpz;
}

template<typename Bigint>
mpz_class bigint_to_mpz(const Bigint &x)
{
    mpz_class x_mpz;

    x.to_mpz(x_mpz.get_mpz_t());

    return x_mpz;
}

template<typename FieldT>
int legendre(const FieldT &x)
{
    static const auto p = bigint_to_mpz(FieldT::mod);

    return mpz_legendre(field_to_mpz(x).get_mpz_t(), p.get_mpz_t());
}

template<typename FieldT>
std::pair<FieldT, FieldT> get_irreducible_pair()
{
    FieldT a{1};

    for (size_t i = 0; i < 4096; ++i, ++a)
    {
        FieldT a2{a.squared()};
        FieldT b{1};

        for (size_t j = 0; j < 4096; ++j, ++b)
            if (legendre(a2 - (b + b + b + b)) < 0)
                return {a, b};
    }

    return {0, 0};
}

template<typename FieldT>
FieldT modular_inverse(const FieldT &x, const FieldT &modulus)
{
    static constexpr mp_size_t n = FieldT::num_limbs;

    using Bigint = libff::bigint<n>;

    if (x.is_zero())
        return 0;

    Bigint mod{modulus.is_zero() ? FieldT::mod : modulus.as_bigint()};
    Bigint u{x.as_bigint()};
    Bigint v{mod};      // both source operands are destroyed by mpn_gcdext
    mp_limb_t g[n];     /* gp should have room for vn = n limbs */
    mp_limb_t s[n + 1]; /* sp should have room for vn+1 limbs */
    mp_size_t sn;

    /* computes gcd(u, v) = g = u*s + v*t, so u*s will be 1 (mod v) */
    mp_size_t gn = mpn_gcdext(g, s, &sn, u.data, n, v.data, n);
    mp_size_t asn = std::abs(sn);

    if (gn != 1 || g[0] != 1) // no inverse
        return 0;

    if (asn >= n)
    {
        /* if sn could require modulus reduction, do it here */
        mpn_tdiv_qr(g, u.data, 0, s, asn, mod.data, n);
    }
    else
    {
        /* otherwise just copy it over */
        mpn_zero(u.data, n);
        mpn_copyi(u.data, s, asn);
    }

    /* fix up the negative sn */
    if (sn < 0 && mpn_sub_n(u.data, mod.data, u.data, n) != 0)
        return 0;

    return u;
}

template<typename FieldT, size_t sz>
static constexpr std::array<FieldT, sz> random_array()
{
    std::array<FieldT, sz> a;

    std::generate(a.begin(), a.end(), []() { return FieldT::random_element(); });

    return a;
}

namespace libff
{
    template<mp_size_t n, const bigint<n> &m>
    Fp_model<n, m> operator++(Fp_model<n, m> &x)
    {
        return x += 1;
    }

    template<mp_size_t n, const bigint<n> &m>
    Fp_model<n, m> operator++(Fp_model<n, m> &x, int)
    {
        auto t{x};
        ++x;

        return t;
    }

    template<mp_size_t n, const bigint<n> &m>
    Fp_model<n, m> operator--(Fp_model<n, m> &x)
    {
        return x -= 1;
    }

    template<mp_size_t n, const bigint<n> &m>
    Fp_model<n, m> operator--(Fp_model<n, m> &x, int)
    {
        auto t{x};
        --x;

        return t;
    }
}; // namespace libff
