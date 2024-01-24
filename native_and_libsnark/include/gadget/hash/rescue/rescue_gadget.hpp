#pragma once

#include "gadget/field_variable.hpp"
#include "gadget/pb_variable_pp.hpp"

template<typename Rescue>
class RescueGadget : public GadgetPP<typename Rescue::Field>
{
public:
    using Field = typename Rescue::Field;
    using Hash = Rescue;
    using DigVar = FieldVariable<Field>;
    using BlockVar = std::array<DigVar, Hash::RATE>;

    static constexpr size_t RATE = Hash::RATE;
    static constexpr size_t CAPACITY = Hash::CAPACITY;
    static constexpr size_t BRANCH_N = Hash::BRANCH_N;
    static constexpr size_t ROUNDS_N = Hash::ROUNDS_N;
    static constexpr size_t DIGEST_SIZE = Hash::DIGEST_SIZE;
    static constexpr size_t BLOCK_SIZE = Hash::BLOCK_SIZE;
    static constexpr size_t DIGEST_VARS = 1;

    const BlockVar in;
    const DigVar out;

private:
    using super = GadgetPP<typename Rescue::Field>;
    using LC = typename super::LC;

    using super::constrain;
    using super::val;

    static constexpr size_t INTER_N = 6 * ROUNDS_N;

    static constexpr auto &alpha = Hash::alpha;
    static constexpr auto &alpha_i = Hash::alpha_i;
    static constexpr auto &mat = Hash::mat;
    static constexpr auto &rc = Hash::round_c;

    std::vector<PbVariablePP<Field>> inter[BRANCH_N];

public:
    static size_t size() { return INTER_N * BRANCH_N; }

    RescueGadget(libsnark::protoboard<Field> &pb, const BlockVar &in, const DigVar &out,
                 const std::string &annotation_prefix) :
        super{pb, annotation_prefix}, in{in}, out{out}
    {
        for (size_t i = 0; i < BRANCH_N; ++i)
            for (size_t j = 0; j < INTER_N; ++j)
                inter[i].emplace_back(pb, FMT(""));
    }

    void generate_r1cs_constraints()
    {
        size_t i[BRANCH_N]{};
        std::array<LC, BRANCH_N> t{};

        for (size_t j = 0; j < RATE; ++j)
            t[j] = in[j][0];

        for (size_t j = 0; j < ROUNDS_N; ++j)
        {
            // Direct SBOX (x -> x^a), assume a == 5
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                // x^2
                i[k] += constrain(t[k], t[k], inter[k][i[k]]);
                // x^4
                i[k] += constrain(inter[k][i[k] - 1], inter[k][i[k] - 1], inter[k][i[k]]);
                // x^5
                i[k] += constrain(inter[k][i[k] - 1], t[k], inter[k][i[k]]);
            }

            // reset temporaries
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = 0;

            // MDS MULTIPLICATION
            for (size_t k = 0; k < BRANCH_N; ++k)
                for (size_t l = 0; l < BRANCH_N; ++l)
                    t[k] = t[k] + mat[k * BRANCH_N + l] * inter[l][i[l] - 1];

            // ADD FIRST CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = t[k] + rc[j * 2 * BRANCH_N + k];

            // Inverse SBOX (x -> x^{1/a}), assume a == 5
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                // x^{1/a} == y    <==>    y^a == x
                // y^2
                ++i[k];
                i[k] += constrain(inter[k][i[k] - 1], inter[k][i[k] - 1], inter[k][i[k]]);
                // y^4
                i[k] += constrain(inter[k][i[k] - 1], inter[k][i[k] - 1], inter[k][i[k]]);
                // y^5
                constrain(inter[k][i[k] - 1], inter[k][i[k] - 3], t[k]);
            }

            // reset temporaries
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = 0;

            // MDS MULTIPLICATION
            for (size_t k = 0; k < BRANCH_N; ++k)
                for (size_t l = 0; l < BRANCH_N; ++l)
                    t[k] = t[k] + mat[k * BRANCH_N + l] * inter[l][i[l] - 3];

            // ADD SECOND CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = t[k] + rc[BRANCH_N + j * 2 * BRANCH_N + k];
        }

        // bind output
        constrain(t[0], 1, out[0]);
    }

    void generate_r1cs_witness()
    {
        size_t i[BRANCH_N]{};
        std::array<Field, BRANCH_N> t{};

        for (size_t j = 0; j < RATE; ++j)
            t[j] = val(in[j][0]);

        for (size_t j = 0; j < ROUNDS_N; ++j)
        {
            // Direct SBOX (x -> x^a), assume a == 5
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                // x^2
                val(inter[k][i[k]]) = t[k] * t[k];
                ++i[k];
                // x^4
                val(inter[k][i[k]]) = val(inter[k][i[k] - 1]) * val(inter[k][i[k] - 1]);
                ++i[k];
                // x^5
                val(inter[k][i[k]]) = val(inter[k][i[k] - 1]) * t[k];
                ++i[k];
            }

            // reset temporaries
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = 0;

            // MDS MULTIPLICATION
            for (size_t k = 0; k < BRANCH_N; ++k)
                for (size_t l = 0; l < BRANCH_N; ++l)
                    t[k] = t[k] + mat[k * BRANCH_N + l] * val(inter[l][i[l] - 1]);

            // ADD FIRST CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = t[k] + rc[j * 2 * BRANCH_N + k];

            // Inverse SBOX (x -> x^{1/a}), assume a == 5
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                // x^{1/a} == y    <==>    y^a == x
                val(inter[k][i[k]]) = t[k];
                Hash::raise_alpha_inv(val(inter[k][i[k]]));
                ++i[k];
                // y^2
                val(inter[k][i[k]]) = val(inter[k][i[k] - 1]) * val(inter[k][i[k] - 1]);
                ++i[k];
                // y^4
                val(inter[k][i[k]]) = val(inter[k][i[k] - 1]) * val(inter[k][i[k] - 1]);
                ++i[k];
            }

            // reset temporaries
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = 0;

            // MDS MULTIPLICATION
            for (size_t k = 0; k < BRANCH_N; ++k)
                for (size_t l = 0; l < BRANCH_N; ++l)
                    t[k] = t[k] + mat[k * BRANCH_N + l] * val(inter[l][i[l] - 3]);

            // ADD SECOND CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = t[k] + rc[BRANCH_N + j * 2 * BRANCH_N + k];
        }

        // bind output
        val(out[0]) = t[0];
    }
};
