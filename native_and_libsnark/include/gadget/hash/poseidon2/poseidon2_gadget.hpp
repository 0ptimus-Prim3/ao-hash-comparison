#pragma once

#include "gadget/field_variable.hpp"
#include "gadget/pb_variable_pp.hpp"

template<typename Poseidon2>
class Poseidon2Gadget : public GadgetPP<typename Poseidon2::Field>
{
public:
    using Field = typename Poseidon2::Field;
    using Hash = Poseidon2;
    using DigVar = FieldVariable<Field>;
    using BlockVar = std::array<DigVar, Hash::BRANCH_N>;

    static constexpr size_t BRANCH_N = Hash::BRANCH_N;
    static constexpr size_t ROUNDS_f = Hash::ROUNDS_f_N;
    static constexpr size_t ROUNDS_F = Hash::ROUNDS_F_N;
    static constexpr size_t ROUNDS_P = Hash::ROUNDS_P_N;
    static constexpr size_t ROUNDS_N = Hash::ROUNDS_N;
    static constexpr size_t DIGEST_SIZE = Hash::DIGEST_SIZE;
    static constexpr size_t BLOCK_SIZE = Hash::BLOCK_SIZE;
    static constexpr size_t DIGEST_VARS = 1;

    const BlockVar in;
    const DigVar out;

private:
    using super = GadgetPP<typename Poseidon2::Field>;
    using LC = typename super::LC;

    using super::constrain;
    using super::val;

    static constexpr size_t INTER0_N = 3 * ROUNDS_N;
    static constexpr size_t INTERk_N = 3 * ROUNDS_F;

    static constexpr auto &ext_rc = Hash::ext_round_c;
    static constexpr auto &int_rc = Hash::int_round_c;
    static constexpr auto &int_mat = Hash::int_mat;

    std::vector<PbVariablePP<Field>> inter[BRANCH_N];

public:
    static size_t size() { return INTER0_N + (BRANCH_N - 1) * INTERk_N; }

    Poseidon2Gadget(libsnark::protoboard<Field> &pb, const BlockVar &in, const DigVar &out,
                    const std::string &annotation_prefix) :
        super{pb, annotation_prefix},
        in{in}, out{out}
    {
        for (size_t i = 0; i < INTER0_N; ++i)
            inter[0].emplace_back(pb, FMT(""));

        for (size_t i = 1; i < BRANCH_N; ++i)
            for (size_t j = 0; j < INTERk_N; ++j)
                inter[i].emplace_back(pb, FMT(""));
    }

#define EXT_MDS_MUL()                                                                              \
    do                                                                                             \
    {                                                                                              \
        if constexpr (BRANCH_N == 1)                                                               \
        {                                                                                          \
            ;                                                                                      \
        }                                                                                          \
        else if constexpr (BRANCH_N == 2)                                                          \
        {                                                                                          \
            s[0] = t[0] + t[1];                                                                    \
            t[0] = t[0] + s[0];                                                                    \
            t[1] = t[1] + t[1];                                                                    \
            t[1] = t[1] + s[0];                                                                    \
        }                                                                                          \
        else if constexpr (BRANCH_N == 3)                                                          \
        {                                                                                          \
            s[0] = t[0] + t[1] + t[2];                                                             \
                                                                                                   \
            t[0] = t[0] + s[0];                                                                    \
            t[1] = t[1] + t[1];                                                                    \
            t[1] = t[1] + s[0];                                                                    \
            t[2] = t[2] + t[2];                                                                    \
            t[2] = t[2] + t[2];                                                                    \
            t[2] = t[2] + s[0];                                                                    \
        }                                                                                          \
        else                                                                                       \
        {                                                                                          \
            for (size_t j = 0; j < BRANCH_N; j += 4)                                               \
            {                                                                                      \
                s[0] = t[j + 0] + t[j + 1];                                                        \
                s[1] = t[j + 2] + t[j + 3];                                                        \
                s[2] = t[j + 1] + t[j + 1] + s[1];                                                 \
                s[3] = t[j + 3] + t[j + 3] + s[0];                                                 \
                                                                                                   \
                t[j + 3] = s[1] + s[1] + s[1] + s[1] + s[3];                                       \
                t[j + 1] = s[0] + s[0] + s[0] + s[0] + s[2];                                       \
                t[j + 0] = s[3] + t[j + 1];                                                        \
                t[j + 2] = s[2] + t[j + 3];                                                        \
            }                                                                                      \
                                                                                                   \
            for (size_t j = 0; j < BRANCH_N; ++j)                                                  \
            {                                                                                      \
                s[j] = t[j] + t[j];                                                                \
                                                                                                   \
                for (size_t k = j & 3; k < j; k += 4)                                              \
                    s[j] = s[j] + t[k];                                                            \
                                                                                                   \
                for (size_t k = j + 4; k < BRANCH_N; k += 4)                                       \
                    s[j] = s[j] + t[k];                                                            \
            }                                                                                      \
                                                                                                   \
            t = s;                                                                                 \
        }                                                                                          \
    } while (0)

#define INT_MDS_MUL()                                                                              \
    do                                                                                             \
    {                                                                                              \
        if constexpr (BRANCH_N <= 3)                                                               \
            EXT_MDS_MUL();                                                                         \
        else                                                                                       \
        {                                                                                          \
            s[0] = t[0];                                                                           \
                                                                                                   \
            for (size_t j = 1; j < BRANCH_N; ++j)                                                  \
                s[0] = s[0] + t[j];                                                                \
                                                                                                   \
            for (size_t j = 0; j < BRANCH_N; ++j)                                                  \
            {                                                                                      \
                t[j] = t[j] * int_mat[j];                                                          \
                t[j] = t[j] + s[0];                                                                \
            }                                                                                      \
        }                                                                                          \
    } while (0)

    void generate_r1cs_constraints()
    {
        std::array<LC, BRANCH_N> s{};
        std::array<LC, BRANCH_N> t{};
        std::array<size_t, BRANCH_N> i{};

        for (size_t j = 0; j < BRANCH_N; ++j)
            t[j] = in[j][0];

        // first matrix multiplication
        EXT_MDS_MUL();

        // INITIAL FULL LAYERS
        for (size_t j = 0; j < ROUNDS_f; ++j)
        {
            // ADD CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = t[k] + ext_rc[j * BRANCH_N + k];

            // FULL SBOX (x^5)
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                i[k] += constrain(t[k], t[k], inter[k][i[k]]);
                i[k] += constrain(inter[k][i[k] - 1], inter[k][i[k] - 1], inter[k][i[k]]);
                i[k] += constrain(t[k], inter[k][i[k] - 1], inter[k][i[k]]);
            }

            // MDS multiplication
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = inter[k][i[k] - 1];
            EXT_MDS_MUL();
        }

        // PARTIAL LAYERS
        for (size_t j = 0; j < ROUNDS_P; ++j)
        {
            // ADD CONSTANTS
            t[0] = t[0] + int_rc[j];

            // SBOX (x^5) for first element
            i[0] += constrain(t[0], t[0], inter[0][i[0]]);
            i[0] += constrain(inter[0][i[0] - 1], inter[0][i[0] - 1], inter[0][i[0]]);
            i[0] += constrain(inter[0][i[0] - 1], t[0], inter[0][i[0]]);

            // MDS multiplication
            t[0] = inter[0][i[0] - 1];
            INT_MDS_MUL();
        }

        // FINAL FULL LAYERS
        for (size_t j = 0; j < ROUNDS_f; ++j)
        {
            // ADD CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = t[k] + ext_rc[(ROUNDS_f + j) * BRANCH_N + k];

            // FULL SBOX (x^5)
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                i[k] += constrain(t[k], t[k], inter[k][i[k]]);
                i[k] += constrain(inter[k][i[k] - 1], inter[k][i[k] - 1], inter[k][i[k]]);
                i[k] += constrain(t[k], inter[k][i[k] - 1], inter[k][i[k]]);
            }

            // MDS multiplication
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = inter[k][i[k] - 1];
            EXT_MDS_MUL();
        }

        constrain(t[0] + in[0][0], 1, out[0]);
    }

    void generate_r1cs_witness()
    {
        std::array<Field, BRANCH_N> s{};
        std::array<Field, BRANCH_N> t{};
        std::array<size_t, BRANCH_N> i{};

        for (size_t j = 0; j < BRANCH_N; ++j)
            t[j] = val(in[j][0]);

        // first matrix multiplication
        EXT_MDS_MUL();

        // INITIAL FULL LAYERS
        for (size_t j = 0; j < ROUNDS_f; ++j)
        {
            // ADD CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] += ext_rc[j * BRANCH_N + k];

            // FULL SBOX (x^5)
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                val(inter[k][i[k]]) = t[k] * t[k];
                ++i[k];
                val(inter[k][i[k]]) = val(inter[k][i[k] - 1]) * val(inter[k][i[k] - 1]);
                ++i[k];
                val(inter[k][i[k]]) = t[k] * val(inter[k][i[k] - 1]);
                ++i[k];
            }

            // MDS multiplication
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = val(inter[k][i[k] - 1]);
            EXT_MDS_MUL();
        }

        // PARTIAL LAYERS
        for (size_t j = 0; j < ROUNDS_P; ++j)
        {
            // ADD CONSTANTS
            t[0] += int_rc[j];

            // SBOX (x^5) for first element
            val(inter[0][i[0]]) = t[0] * t[0];
            ++i[0];
            val(inter[0][i[0]]) = val(inter[0][i[0] - 1]) * val(inter[0][i[0] - 1]);
            ++i[0];
            val(inter[0][i[0]]) = t[0] * val(inter[0][i[0] - 1]);
            ++i[0];

            // MDS multiplication
            t[0] = val(inter[0][i[0] - 1]);
            INT_MDS_MUL();
        }

        // FINAL FULL LAYERS
        for (size_t j = 0; j < ROUNDS_f; ++j)
        {
            // ADD CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] += ext_rc[(ROUNDS_f + j) * BRANCH_N + k];

            // FULL SBOX (x^5)
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                val(inter[k][i[k]]) = t[k] * t[k];
                ++i[k];
                val(inter[k][i[k]]) = val(inter[k][i[k] - 1]) * val(inter[k][i[k] - 1]);
                ++i[k];
                val(inter[k][i[k]]) = t[k] * val(inter[k][i[k] - 1]);
                ++i[k];
            }

            // MDS multiplication
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = val(inter[k][i[k] - 1]);
            EXT_MDS_MUL();
        }

        val(out[0]) = t[0] + val(in[0][0]);
    }

#undef INT_MDS_MUL
#undef EXT_MDS_MUL
};
