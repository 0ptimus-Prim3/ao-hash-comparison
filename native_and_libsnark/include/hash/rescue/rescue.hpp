#pragma once

#include "util/algebra.hpp"

template<typename FieldT = libff::Fr<libff::default_ec_pp>, size_t rate = 2, size_t capacity = 2,
         size_t rounds = 19>
class Rescue
{
public:
    using Field = FieldT;

    static constexpr size_t RATE = rate;
    static constexpr size_t CAPACITY = capacity;
    static constexpr size_t ROUNDS_N = rounds;
    static constexpr size_t DIGEST_SIZE = field_size<Field>();
    static constexpr size_t BLOCK_SIZE = DIGEST_SIZE * RATE;
    static constexpr size_t BRANCH_N = RATE + CAPACITY;

    using Sponge = std::array<Field, BRANCH_N>;
    using Constants = std::array<Field, ROUNDS_N * 2 * BRANCH_N>;
    using Matrix = std::array<Field, BRANCH_N * BRANCH_N>;

    static inline const struct Init
    {
        Init() { libff::default_ec_pp::init_public_params(); }
    } init;

    static inline const Field alpha{5};
    static inline const Field alpha_i{modular_inverse(alpha, Field{-1})};

    static inline Constants gen_constants()
    {
        Constants c{};

        std::iota(c.begin(), c.end(), Field{1});

        return c;
    }

    static inline Matrix gen_matrix()
    {
        Matrix m;

        std::iota(m.begin(), m.end(), Field{1});

        return m;
    }

    static inline const Constants round_c{gen_constants()};
    static inline const Matrix mat{gen_matrix()};

    static void raise_alpha(Field &x)
    {
        Field t{x};

        x *= x;
        x *= x;
        x *= t;
    }

    static void raise_alpha_inv(Field &x)
    {
        static const auto ai = alpha_i.as_bigint();

        x ^= ai;
    }

    static void matmul(Sponge &x)
    {
        Sponge s{};
        for (size_t i = 0; i < BRANCH_N; ++i)
            for (size_t j = 0; j < BRANCH_N; ++j)
                s[i] += mat[i * BRANCH_N + j] * x[j];

        x = s;
    }

    static void hash_field(Sponge &h)
    {
        for (size_t i = 0; i < ROUNDS_N; ++i)
        {
            // Direct SBOX
            for (size_t j = 0; j < BRANCH_N; ++j)
                raise_alpha(h[j]);

            // MDS
            matmul(h);

            // ADD FIRST CONSTANTS
            for (size_t j = 0; j < BRANCH_N; ++j)
                h[j] += round_c[i * 2 * BRANCH_N + j];

            // Inverse SBOX
            for (size_t j = 0; j < BRANCH_N; ++j)
                raise_alpha_inv(h[j]);

            // Second MDS
            matmul(h);

            // ADD SECOND CONSTANTS
            for (size_t j = 0; j < BRANCH_N; ++j)
                h[j] += round_c[BRANCH_N + i * 2 * BRANCH_N + j];
        }
    }

    static void hash_oneblock(uint8_t *digest, const void *message)
    {
        Sponge h{};
        mpz_class tmp;

        for (size_t i = 0; i < RATE; ++i)
        {
            mpz_import(tmp.get_mpz_t(), DIGEST_SIZE, 1, 1, 0, 0,
                       (const char *)message + DIGEST_SIZE * i);

            h[i] = Field{tmp.get_mpz_t()};
        }

        hash_field(h);
        h[0].as_bigint().to_mpz(tmp.get_mpz_t());

        memset(digest, 0, DIGEST_SIZE);
        mpz_export(digest, NULL, 1, 1, 0, 0, tmp.get_mpz_t());
    }

    static void hash_add(void *x, const void *y)
    {
        mpz_class tmp;

        mpz_import(tmp.get_mpz_t(), DIGEST_SIZE, 1, 1, 0, 0, x);
        Field xf{tmp.get_mpz_t()};

        mpz_import(tmp.get_mpz_t(), DIGEST_SIZE, 1, 1, 0, 0, y);
        Field yf{tmp.get_mpz_t()};

        xf += yf;

        xf.as_bigint().to_mpz(tmp.get_mpz_t());
        memset(x, 0, DIGEST_SIZE);
        mpz_export(x, NULL, 1, 1, 0, 0, tmp.get_mpz_t());
    }

    Rescue() = delete;
};

// Valid realizations of a template class must be initialized before main()!
template class Rescue<>;
