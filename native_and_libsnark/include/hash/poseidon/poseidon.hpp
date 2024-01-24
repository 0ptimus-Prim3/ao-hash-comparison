#pragma once

#include "util/algebra.hpp"

template<typename FieldT = libff::Fr<libff::default_ec_pp>, size_t rate = 2, size_t capacity = 1,
         size_t rounds_f = 4, size_t rounds_p = 57>
class Poseidon
{
public:
    using Field = FieldT;

    static constexpr size_t RATE = rate;
    static constexpr size_t CAPACITY = capacity;
    static constexpr size_t ROUNDS_f_N = rounds_f;
    static constexpr size_t ROUNDS_P_N = rounds_p;
    static constexpr size_t DIGEST_SIZE = field_size<Field>();
    static constexpr size_t BLOCK_SIZE = DIGEST_SIZE * RATE;
    static constexpr size_t BRANCH_N = RATE + CAPACITY;
    static constexpr size_t ROUNDS_F_N = 2 * ROUNDS_f_N;
    static constexpr size_t ROUNDS_N = ROUNDS_F_N + ROUNDS_P_N;
    static constexpr size_t CONST_N = BRANCH_N * ROUNDS_N;

    using Sponge = std::array<Field, BRANCH_N>;
    using Matrix = std::array<Field, BRANCH_N * BRANCH_N>;
    using Constants = std::array<Field, ROUNDS_N * BRANCH_N>;

    static inline const struct Init
    {
        Init() { field_init<Field>(); }
    } init;

    static inline const Constants round_c{[]
                                          {
                                              Constants c;

                                              std::iota(c.begin(), c.end(), Field{1});

                                              return c;
                                          }()};
    static inline const Matrix mds_mat{[]
                                       {
                                           Matrix m;
                                           std::array<Field, BRANCH_N> x;
                                           std::array<Field, BRANCH_N> y;

                                           std::iota(x.begin(), x.end(), Field{1});
                                           std::iota(y.begin(), y.end(), Field{BRANCH_N + 1});

                                           for (size_t i = 0; i < BRANCH_N; ++i)
                                               for (size_t j = 0; j < BRANCH_N; ++j)
                                                   m[i * BRANCH_N + j] = field_inverse(x[i] + y[j]);

                                           return m;
                                       }()};

    static void fifth(Field &x)
    {
        Field t{x};

        x *= x;
        x *= x;
        x *= t;
    }

    static void matmul(Sponge &arr)
    {
        Sponge sum{};

        for (size_t i = 0; i < BRANCH_N; ++i)
            for (size_t j = 0; j < BRANCH_N; ++j)
                sum[i] += mds_mat[i * BRANCH_N + j] * arr[j];

        arr = sum;
    }

    static Field hash_field(Sponge &h)
    {
        // INITIAL FULL LAYERS
        for (size_t i = 0; i < ROUNDS_f_N; ++i)
        {
            for (size_t j = 0; j < BRANCH_N; ++j)
                h[j] += round_c[i * BRANCH_N + j];

            for (size_t j = 0; j < BRANCH_N; ++j)
                fifth(h[j]);

            matmul(h);
        }


        // PARTIAL LAYERS
        for (size_t i = 0; i < ROUNDS_P_N; ++i)
        {
            for (size_t j = 0; j < BRANCH_N; ++j)
                h[j] += round_c[(ROUNDS_f_N + i) * BRANCH_N + j];

            fifth(h[0]);
            matmul(h);
        }

        // FINAL FULL LAYERS
        for (size_t i = 0; i < ROUNDS_f_N; ++i)
        {
            for (size_t j = 0; j < BRANCH_N; ++j)
                h[j] += round_c[(ROUNDS_f_N + ROUNDS_P_N + i) * BRANCH_N + j];

            for (size_t j = 0; j < BRANCH_N; ++j)
                fifth(h[j]);

            matmul(h);
        }

        return h[0];
    }

    static void hash_oneblock(uint8_t *digest, const void *message)
    {
        Sponge h{};

        field_load(h.data(), message, RATE);
        hash_field(h);
        field_store(digest, h.data(), 1);
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

    Poseidon() = delete;
};

// Valid realizations of a template class must be initialized before main()!
template class Poseidon<>;