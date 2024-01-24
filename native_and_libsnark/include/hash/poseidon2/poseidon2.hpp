#pragma once

#include "util/algebra.hpp"
#include "util/string_utils.hpp"


template<typename FieldT = libff::Fr<libff::default_ec_pp>, size_t branch_n = 2,
         size_t rounds_f = 4, size_t rounds_p = 57>
class Poseidon2
{
public:
    static_assert(branch_n < 4 || branch_n % 4 == 0, "Invalid branch size");
    using Field = FieldT;

    static constexpr size_t BRANCH_N = branch_n;
    static constexpr size_t ROUNDS_f_N = rounds_f;
    static constexpr size_t ROUNDS_P_N = rounds_p;
    static constexpr size_t FIELD_SIZE = field_size<Field>();
    static constexpr size_t DIGEST_SIZE = FIELD_SIZE;
    static constexpr size_t BLOCK_SIZE = DIGEST_SIZE * BRANCH_N;
    static constexpr size_t ROUNDS_F_N = 2 * ROUNDS_f_N;
    static constexpr size_t ROUNDS_N = ROUNDS_F_N + ROUNDS_P_N;
    static constexpr size_t EXT_CONST_N = BRANCH_N * ROUNDS_F_N;

    using Block = std::array<Field, BRANCH_N>;
    using IntMatrix = std::array<Field, BRANCH_N>;
    using ExtConstants = std::array<Field, EXT_CONST_N>;
    using IntConstants = std::array<Field, ROUNDS_P_N>;

    static inline const struct Init
    {
        Init() { field_init<Field>(); }
    } init;

    static inline const ExtConstants ext_round_c{[]
                                                 {
                                                     ExtConstants c;

                                                     std::iota(c.begin(), c.end(), Field{1});

                                                     return c;
                                                 }()};

    static inline const IntConstants int_round_c{[]
                                                 {
                                                     IntConstants c;

                                                     std::iota(c.begin(), c.end(),
                                                               Field{EXT_CONST_N + 1});

                                                     return c;
                                                 }()};

    static inline const IntMatrix int_mat{[]
                                          {
                                              IntMatrix m;

                                              std::iota(m.begin(), m.end(), Field{1});

                                              return m;
                                          }()};

    static void fifth(Field &x)
    {
        Field t{x};

        x *= x;
        x *= x;
        x *= t;
    }

    static void ext_matmul(Block &x)
    {
        if constexpr (BRANCH_N == 1)
            return;

        if constexpr (BRANCH_N == 2)
        {
            // M_E = [2, 1; 3, 1]
            Field s = x[0] + x[1];

            x[0] += s;
            x[1] += x[1];
            x[1] += s;
        }
        else if constexpr (BRANCH_N == 3)
        {
            // M_E = [2, 1, 1; 1, 3, 1; 1, 1, 5]
            Field s = x[0] + x[1] + x[2];

            x[0] += s;
            x[1] += x[1];
            x[1] += s;
            x[2] += x[2];
            x[2] += x[2];
            x[2] += s;
        }
        else
        {
            Block t;

            for (size_t i = 0; i < BRANCH_N; i += 4)
            {
                t[0] = x[i + 0] + x[i + 1];
                t[1] = x[i + 2] + x[i + 3];
                t[2] = x[i + 1] + x[i + 1] + t[1];
                t[3] = x[i + 3] + x[i + 3] + t[0];

                x[i + 3] = t[1] + t[1] + t[1] + t[1] + t[3];
                x[i + 1] = t[0] + t[0] + t[0] + t[0] + t[2];
                x[i + 0] = t[3] + x[i + 1];
                x[i + 2] = t[2] + x[i + 3];
            }

            for (size_t i = 0; i < BRANCH_N; ++i)
            {
                t[i] = x[i] + x[i];

                for (size_t j = i & 3; j < i; j += 4)
                    t[i] += x[j];

                for (size_t j = i + 4; j < BRANCH_N; j += 4)
                    t[i] += x[j];
            }

            x = t;
        }
    }

    static void int_matmul(Block &x)
    {
        if constexpr (BRANCH_N <= 3)
            ext_matmul(x);
        else
        {
            Field s{0ULL};

            for (size_t i = 0; i < BRANCH_N; ++i)
                s += x[i];

            for (size_t i = 0; i < BRANCH_N; ++i)
            {
                x[i] *= int_mat[i];
                x[i] += s;
            }
        }
    }


    static void hash_field(Block &x)
    {
        Field t{x[0]};

        // first matrix multiplication
        ext_matmul(x);

        // INITIAL FULL LAYERS
        for (size_t i = 0; i < ROUNDS_f_N; ++i)
        {
            for (size_t j = 0; j < BRANCH_N; ++j)
                x[j] += ext_round_c[i * BRANCH_N + j];

            for (size_t j = 0; j < BRANCH_N; ++j)
                fifth(x[j]);

            ext_matmul(x);
        }

        // PARTIAL LAYERS
        for (size_t i = 0; i < ROUNDS_P_N; ++i)
        {
            x[0] += int_round_c[i];

            fifth(x[0]);
            int_matmul(x);
        }

        // FINAL FULL LAYERS
        for (size_t i = 0; i < ROUNDS_f_N; ++i)
        {
            for (size_t j = 0; j < BRANCH_N; ++j)
                x[j] += ext_round_c[(ROUNDS_f_N + i) * BRANCH_N + j];

            for (size_t j = 0; j < BRANCH_N; ++j)
                fifth(x[j]);

            ext_matmul(x);
        }

        x[0] += t;
    }

    static void hash_oneblock(uint8_t *digest, const void *message)
    {
        Block x{};

        field_load(x.data(), message, BRANCH_N);
        hash_field(x);
        field_store(digest, x.data(), 1);
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

    Poseidon2() = delete;
};

// Valid realizations of a template class must be initialized before main()!
template class Poseidon2<>;
