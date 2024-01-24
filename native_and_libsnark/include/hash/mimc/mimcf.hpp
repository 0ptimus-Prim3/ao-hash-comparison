#pragma once

#include "util/algebra.hpp"

template<typename FieldT = libff::Fr<libff::default_ec_pp>, size_t rounds = 320>
class MimcF
{
public:
    using Field = FieldT;

    static constexpr size_t ROUNDS_N = rounds;
    static constexpr size_t DIGEST_SIZE = field_size<Field>();
    static constexpr size_t BLOCK_SIZE = 2 * DIGEST_SIZE;

    using Feistel = std::array<Field, 2>;
    using Constants = std::array<Field, ROUNDS_N>;

    static inline const struct Init
    {
        Init() { libff::default_ec_pp::init_public_params(); }
    } init;

    static Constants gen_roundc()
    {
        Constants c;

        std::iota(c.begin(), c.end(), Field{0});

        return c;
    }

    static inline const Constants round_c{gen_roundc()};
    static inline const Field d{3};

    static void cube(FieldT &x)
    {
        FieldT t{x};

        x *= x;
        x *= t;
    }

    static void hash_field(Feistel &x) {
        
    }

    static void hash_oneblock(uint8_t *digest, const void *message)
    {
        Feistel h{};
        mpz_class tmp;

        for (size_t i = 0; i < 2; ++i)
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
        FieldT xf{tmp.get_mpz_t()};

        mpz_import(tmp.get_mpz_t(), DIGEST_SIZE, 1, 1, 0, 0, y);
        FieldT yf{tmp.get_mpz_t()};

        xf += yf;

        xf.as_bigint().to_mpz(tmp.get_mpz_t());
        memset(x, 0, DIGEST_SIZE);
        mpz_export(x, NULL, 1, 1, 0, 0, tmp.get_mpz_t());
    }

    MimcF() = delete;
};

// Valid realizations of a template class must be initialized before main()!
template class MimcF<>;
