#pragma once

#include "util/algebra.hpp"

template<typename FieldT = libff::Fr<libff::default_ec_pp>, size_t rounds = 160>
class Mimc256
{
public:
    using Field = FieldT;
    static constexpr size_t DIGEST_SIZE = field_size<Field>();
    static constexpr size_t BLOCK_SIZE = 2 * DIGEST_SIZE;
    static constexpr size_t ROUNDS_N = rounds;

    static inline const struct Init
    {
        Init() { libff::default_ec_pp::init_public_params(); }
    } init;

    static std::array<FieldT, ROUNDS_N - 1> gen_roundc()
    {
        std::array<FieldT, ROUNDS_N - 1> c;

        std::iota(c.begin(), c.end(), Field{1});

        return c;
    }


    static inline const std::array<FieldT, ROUNDS_N - 1> round_c{gen_roundc()};

    static void cube(FieldT &x)
    {
        FieldT t{x};
        x *= t;
        x *= t;
    }

    static FieldT hash_field(const FieldT &x, const FieldT &y)
    {
        // optimize first round
        FieldT h{x};

        cube(h);

        for (size_t i = 0; i < ROUNDS_N - 1; ++i)
        {
            h += round_c[i];
            cube(h);
        }

        // We use the Davies-Meyer construction
        // optimize first round
        h += y;
        cube(h);

        for (size_t i = 0; i < ROUNDS_N - 1; ++i)
        {
            h += y;
            h += round_c[i];
            cube(h);
        }
        h += y;

        return h;
    }

    static void hash_oneblock(uint8_t *digest, const void *message)
    {
        mpz_class tmp;

        mpz_import(tmp.get_mpz_t(), DIGEST_SIZE, 1, 1, 0, 0, message);
        FieldT x{tmp.get_mpz_t()};

        mpz_import(tmp.get_mpz_t(), DIGEST_SIZE, 1, 1, 0, 0, (const char *)message + DIGEST_SIZE);
        FieldT y{tmp.get_mpz_t()};

        x = hash_field(x, y);
        x.as_bigint().to_mpz(tmp.get_mpz_t());

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

    Mimc256() = delete;
};

// Valid realizations of a template class must be initialized before main()!
template class Mimc256<>;
