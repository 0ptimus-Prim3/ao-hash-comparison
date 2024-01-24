#pragma once

#include "util/algebra.hpp"

template<typename FieldT = libff::Fr<libff::default_ec_pp>, size_t rate = 2, size_t capacity = 2,
         size_t rounds = 19>
class Anemoi
{
public:
    using Field = FieldT;

    static constexpr size_t RATE = rate;
    static constexpr size_t CAPACITY = capacity;
    static constexpr size_t ROUNDS_N = rounds;
    static constexpr size_t DIGEST_SIZE = field_size<Field>();
    static constexpr size_t BLOCK_SIZE = DIGEST_SIZE * RATE;
    static constexpr size_t BRANCH_N = RATE + CAPACITY;

    static_assert(BRANCH_N % 2 == 0, "Anemoi: branch number (rate+capacity) must be even");

    static constexpr size_t ELL = BRANCH_N / 2;

    using Sponge = std::array<Field, BRANCH_N>;
    using State = std::array<Field, ELL>;
    using Constants = std::array<std::pair<Field, Field>, ROUNDS_N * ELL>;
    using Matrix = std::array<Field, ELL * ELL>;

    static inline const struct Init
    {
        Init() { libff::default_ec_pp::init_public_params(); }
    } init;

    static inline const Field g{7};
    static inline const Field g_i{modular_inverse(g, Field{-1})};
    static inline const Field alpha{5};
    static inline const Field alpha_i{modular_inverse(alpha, Field{-1})};

    static inline Constants gen_constants()
    {
        static constexpr size_t N = std::max(ELL, ROUNDS_N);
        static const Field pi_0{(long)14159265358979323846ULL};
        static const Field pi_1{(long)8214808651328230664ULL};

        Constants c;
        std::array<Field, N> t0;
        std::array<Field, N> t1;
        Field t;

        t0[0] = 1;
        t1[0] = 1;
        for (size_t i = 1; i < N; ++i)
        {
            t0[i] = t0[i - 1] * pi_0;
            t1[i] = t1[i - 1] * pi_1;
        }

        for (size_t i = 0; i < ROUNDS_N; ++i)
        {
            for (size_t j = 0; j < ELL; ++j)
            {
                t = t0[i];
                t += t1[j];
                raise_alpha(t);

                c[i * ELL + j].first = g;
                c[i * ELL + j].first *= t0[i];
                c[i * ELL + j].first *= t0[i];
                c[i * ELL + j].first += t;

                c[i * ELL + j].second = g;
                c[i * ELL + j].second *= t1[j];
                c[i * ELL + j].second *= t1[j];
                c[i * ELL + j].second += t;
                c[i * ELL + j].second += g_i;
            }
        }

        return c;
    }

    static inline Matrix gen_matrix()
    {
        static const Field g2 = g * g;
        static const Field g1 = g + 1;

        if constexpr (ELL == 1)
            return Matrix{1};

        if constexpr (ELL == 2)
            return Matrix{
                1, g,            //
                g, g2 + Field{1} //
            };

        if constexpr (ELL == 3)
            return Matrix{
                g1, 1, g1, //
                1,  1, g,  //
                g,  1, 1   //
            };

        if constexpr (ELL == 4)
            return Matrix{
                1,  g2,     g2, g1,     //
                g1, g + g2, g2, g1 + g, //
                g,  g1,     1,  g,      //
                g,  g1 + g, g1, g1      //
            };

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

    static void matmul(State &x)
    {
        if constexpr (ELL == 1)
        {
            ;
        }
        if constexpr (ELL == 2)
        {
            x[0] += g * x[1];
            x[1] += g * x[0];
        }
        if constexpr (ELL == 3)
        {
            Field t{x[0] + g * x[2]};

            x[2] += x[1];
            x[2] += g * x[0];
            x[0] = t;
            x[0] += x[2];
            x[1] += t;
        }
        else if constexpr (ELL == 4)
        {
            x[0] += x[1];
            x[2] += x[3];
            x[3] += g * x[0];
            x[1] = g * (x[1] + x[2]);
            x[0] += x[1];
            x[2] += g * x[3];
            x[1] += x[2];
            x[3] += x[0];
        }
        else
        {
            State s{};
            for (size_t i = 0; i < ELL; ++i)
                for (size_t j = 0; j < ELL; ++j)
                    s[i] += mat[i * ELL + j] * x[j];

            x = s;
        }
    }

    static void rho(State &x) { std::rotate(x.begin(), x.begin() + 1, x.end()); }

    static void hash_field(Sponge &h)
    {
        State x, y;
        Field t;

        for (size_t i = 0; i < ELL; ++i)
        {
            x[i] = h[i];
            y[i] = h[ELL + i];
        }

        for (size_t i = 0; i < ROUNDS_N; ++i)
        {
            // ADD ROUND CONSTANTS
            for (size_t j = 0; j < ELL; ++j)
            {
                x[j] += round_c[i * ELL + j].first;
                y[j] += round_c[i * ELL + j].second;
            }

            // M_x MULTIPLICATION
            matmul(x);
            // RHO ROTATION
            rho(y);
            // M_y MULTIPLICATION
            matmul(y);

            // PSEUDO-HADAMARD
            for (size_t j = 0; j < ELL; ++j)
            {
                x[j] += y[j];
                y[j] += x[j];
            }

            // Flystel SBOX
            /*
            Flystel performs the following computations:
                1. x_1 = x_0 - (g(y_0 * y_0) + g_i)
                2. y_1 = y_0 - x_1^(a_i)
                3. x_2 = x_1 + g(y_1 * y_1)
            */
            for (size_t j = 0; j < ELL; ++j)
            {
                t = y[j];
                t *= t;
                t *= g;
                t += g_i;
                x[j] -= t;
                t = x[j];
                raise_alpha_inv(t);
                y[j] -= t;
                t = y[j];
                t *= t;
                t *= g;
                x[j] += t;
            }
        }

        // Final matrix multiplication
        matmul(x);
        rho(y);
        matmul(y);

        for (size_t i = 0; i < ELL; ++i)
        {
            h[i] = x[i];
            h[ELL + i] = y[i];
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

    Anemoi() = delete;
};

// Valid realizations of a template class must be initialized before main()!
template class Anemoi<>;
