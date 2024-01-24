#pragma once

#include "gadget/field_variable.hpp"
#include "gadget/pb_variable_pp.hpp"

template<typename Poseidon>
class PoseidonGadget : public GadgetPP<typename Poseidon::Field>
{
public:
    using Field = typename Poseidon::Field;
    using Hash = Poseidon;
    using DigVar = FieldVariable<Field>;
    using BlockVar = std::array<DigVar, Hash::RATE>;

    static constexpr size_t RATE = Hash::RATE;
    static constexpr size_t CAPACITY = Hash::CAPACITY;
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
    using super = GadgetPP<typename Poseidon::Field>;
    using LC = typename super::LC;

    using super::constrain;
    using super::val;

    static constexpr size_t INTER0_N = 3 * ROUNDS_N;
    static constexpr size_t INTERk_N = 3 * ROUNDS_F;

    static constexpr auto &rc = Hash::round_c;
    static constexpr auto &mds = Hash::mds_mat;

    std::vector<PbVariablePP<Field>> inter[BRANCH_N];

public:
    static size_t size() { return INTER0_N + (BRANCH_N - 1) * INTERk_N; }

    PoseidonGadget(libsnark::protoboard<Field> &pb, const BlockVar &in, const DigVar &out,
                   const std::string &annotation_prefix) :
        super{pb, annotation_prefix}, in{in}, out{out}
    {
        for (size_t i = 0; i < INTER0_N; ++i)
            inter[0].emplace_back(pb, FMT(""));

        for (size_t i = 1; i < BRANCH_N; ++i)
            for (size_t j = 0; j < INTERk_N; ++j)
                inter[i].emplace_back(pb, FMT(""));
    }

    void generate_r1cs_constraints()
    {
        LC s[BRANCH_N]{};
        LC t[BRANCH_N]{};
        size_t i[BRANCH_N]{};
        size_t ri = 0;

        for (size_t j = 0; j < RATE; ++j)
            t[j] = in[j][0];

        // INITIAL FULL LAYERS
        for (size_t j = 0; j < ROUNDS_f; ++j)
        {
            // ADD CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = t[k] + rc[ri++];

            // FULL SBOX (x^5)
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                i[k] += constrain(t[k], t[k], inter[k][i[k]]);
                i[k] += constrain(inter[k][i[k] - 1], inter[k][i[k] - 1], inter[k][i[k]]);
                i[k] += constrain(t[k], inter[k][i[k] - 1], inter[k][i[k]]);
            }
            // MDS multiplication
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                s[k] = 0;
                for (size_t l = 0; l < BRANCH_N; ++l)
                    s[k] = s[k] + mds[k * BRANCH_N + l] * inter[l][i[l] - 1];
            }

            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = s[k];
        }

        // PARTIAL LAYERS
        for (size_t j = 0; j < ROUNDS_P; ++j)
        {
            // ADD CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = t[k] + rc[ri++];

            // SBOX (x^5) for first element
            i[0] += constrain(t[0], t[0], inter[0][i[0]]);
            i[0] += constrain(inter[0][i[0] - 1], inter[0][i[0] - 1], inter[0][i[0]]);
            i[0] += constrain(inter[0][i[0] - 1], t[0], inter[0][i[0]]);
            t[0] = inter[0][i[0] - 1];

            // MDS multiplication
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                s[k] = 0;
                for (size_t l = 0; l < BRANCH_N; ++l)
                    s[k] = s[k] + mds[k * BRANCH_N + l] * t[l];
            }

            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = s[k];
        }

        // FINAL FULL LAYERS
        for (size_t j = 0; j < ROUNDS_f; ++j)
        {
            // ADD CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = t[k] + rc[ri++];

            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                i[k] += constrain(t[k], t[k], inter[k][i[k]]);
                i[k] += constrain(inter[k][i[k] - 1], inter[k][i[k] - 1], inter[k][i[k]]);
                i[k] += constrain(t[k], inter[k][i[k] - 1], inter[k][i[k]]);
            }

            // MDS multiplication
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                s[k] = 0;
                for (size_t l = 0; l < BRANCH_N; ++l)
                    s[k] = s[k] + mds[k * BRANCH_N + l] * inter[l][i[l] - 1];
            }

            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = s[k];
        }

        constrain(t[0], 1, out[0]);
    }

    void generate_r1cs_witness()
    {
        Field s[BRANCH_N]{};
        Field t[BRANCH_N]{};
        size_t i[BRANCH_N]{};
        size_t ri = 0;

        for (size_t j = 0; j < RATE; ++j)
            t[j] = val(in[j][0]);

        // INITIAL FULL LAYERS
        for (size_t j = 0; j < ROUNDS_f; ++j)
        {
            // ADD CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] += rc[ri++];

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
            {
                s[k] = 0;
                for (size_t l = 0; l < BRANCH_N; ++l)
                    s[k] += mds[k * BRANCH_N + l] * val(inter[l][i[l] - 1]);
            }

            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = s[k];
        }
        

        // PARTIAL LAYERS
        for (size_t j = 0; j < ROUNDS_P; ++j)
        {
            // ADD CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] += rc[ri++];

            // SBOX (x^5) for first element
            val(inter[0][i[0]]) = t[0] * t[0];
            ++i[0];
            val(inter[0][i[0]]) = val(inter[0][i[0] - 1]) * val(inter[0][i[0] - 1]);
            ++i[0];
            val(inter[0][i[0]]) = t[0] * val(inter[0][i[0] - 1]);
            ++i[0];
            t[0] = val(inter[0][i[0] - 1]);

            // MDS multiplication
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                s[k] = 0;
                for (size_t l = 0; l < BRANCH_N; ++l)
                    s[k] += mds[k * BRANCH_N + l] * t[l];
            }
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = s[k];
        }
        
        // FINAL FULL LAYERS
        for (size_t j = 0; j < ROUNDS_f; ++j)
        {
            // ADD CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] += rc[ri++];

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
            {
                s[k] = 0;
                for (size_t l = 0; l < BRANCH_N; ++l)
                    s[k] += mds[k * BRANCH_N + l] * val(inter[l][i[l] - 1]);
            }

            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = s[k];
        }

        val(out[0]) = t[0];
    }
};
