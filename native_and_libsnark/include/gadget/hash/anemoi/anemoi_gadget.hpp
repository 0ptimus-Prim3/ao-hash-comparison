#pragma once

#include "gadget/field_variable.hpp"
#include "gadget/pb_variable_pp.hpp"

template<typename Anemoi>
class AnemoiGadget : public GadgetPP<typename Anemoi::Field>
{
public:
    using Field = typename Anemoi::Field;
    using Hash = Anemoi;
    using DigVar = FieldVariable<Field>;
    using BlockVar = std::array<DigVar, Hash::RATE>;

    static constexpr size_t RATE = Hash::RATE;
    static constexpr size_t CAPACITY = Hash::CAPACITY;
    static constexpr size_t BRANCH_N = Hash::BRANCH_N;
    static constexpr size_t ELL = Hash::ELL;
    static constexpr size_t ROUNDS_N = Hash::ROUNDS_N;
    static constexpr size_t DIGEST_SIZE = Hash::DIGEST_SIZE;
    static constexpr size_t BLOCK_SIZE = Hash::BLOCK_SIZE;
    static constexpr size_t DIGEST_VARS = 1;

    const BlockVar in;
    const DigVar out;

private:
    using super = GadgetPP<typename Anemoi::Field>;
    using LC = typename super::LC;

    using super::constrain;
    using super::val;

    static constexpr size_t INTER_N = 5 * ROUNDS_N;

    static constexpr auto &alpha = Hash::alpha;
    static constexpr auto &g = Hash::g;
    static constexpr auto &g_i = Hash::g_i;
    static constexpr auto &mat = Hash::mat;
    static constexpr auto &rc = Hash::round_c;

    std::vector<PbVariablePP<Field>> inter[ELL];

public:
    static size_t size() { return ELL * INTER_N; }


    AnemoiGadget(libsnark::protoboard<Field> &pb, const BlockVar &in, const DigVar &out,
                 const std::string &annotation_prefix) :
        super{pb, annotation_prefix}, in{in}, out{out}
    {
        for (size_t i = 0; i < ELL; ++i)
            for (size_t j = 0; j < INTER_N; ++j)
                inter[i].emplace_back(pb, FMT(""));
    }

    void generate_r1cs_constraints()
    {
        size_t i[ELL]{};
        LC tx[ELL]{};
        LC ty[ELL]{};
        LC sx[ELL]{};
        LC sy[ELL]{};

        for (size_t j = 0; j < std::min(RATE, ELL); ++j)
            tx[j] = in[j][0];

        for (size_t j = ELL; j < RATE; ++j)
            ty[j - ELL] = in[j][0];

        for (size_t j = 0; j < ROUNDS_N; ++j)
        {
            // ADD CONSTANTS
            for (size_t k = 0; k < ELL; ++k)
            {
                tx[k] = tx[k] + rc[j * ELL + k].first;
                ty[k] = ty[k] + rc[j * ELL + k].second;
            }

            // M_x MULTIPLICATION
            for (size_t k = 0; k < ELL; ++k)
                for (size_t l = 0; l < ELL; ++l)
                    sx[k] = sx[k] + mat[k * ELL + l] * tx[l];

            // Rho ROTATION
            std::rotate(ty, ty + 1, ty + ELL);

            // M_y MULTIPLICATION
            for (size_t k = 0; k < ELL; ++k)
                for (size_t l = 0; l < ELL; ++l)
                    sy[k] = sy[k] + mat[k * ELL + l] * ty[l];

            // PSEUDO-HADAMARD
            for (size_t k = 0; k < ELL; ++k)
            {
                tx[k] = sx[k] + sy[k];
                ty[k] = sy[k] + tx[k];
            }

            // Flystel SBOX
            /*
            Flystel performs the following computations:
                1. x_1 = x_0 - (g(y_0 * y_0) + g_i)
                2. y_1 = y_0 - x_1^(a_i)
                3. x_2 = x_1 + g(y_1 * y_1)
            */
            for (size_t k = 0; k < ELL; ++k)
            {
                /*
                1.  x_1 = x_0 - (g(y_0 * y_0) + g_i)    <==>
                    x_1 - x_0 + g_i = -gy_0 * y_0
                */
                // tx = x_0, ty = y_0
                // inter[0] = x_1
                i[k] += constrain(-g * ty[k], ty[k], inter[k][i[k]] - tx[k] + g_i);

                // 2. y_1 = y_0 - x_1^(a_i)    <==>    (y_0 - y_1)^a = x_1
                // inter[1] = y_1
                ty[k] = ty[k] - inter[k][i[k]];
                // inter[2] = (y_0 - y_1)^2
                i[k] += constrain(ty[k], ty[k], inter[k][i[k] + 1]);
                // inter[3] = (y_0 - y_1)^4
                i[k] += constrain(inter[k][i[k]], inter[k][i[k]], inter[k][i[k] + 1]);
                // (y_0 - y_1)^5 == inter[0]
                i[k] += constrain(ty[k], inter[k][i[k]], inter[k][i[k] - 3]);

                // 3. x_2 = x_1 + g(y_1 * y_1) <==>    x_2 - x_1 = gy_1 * y_1
                // inter[4] = x_2
                i[k] += constrain(g * inter[k][i[k] - 3], inter[k][i[k] - 3],
                                  inter[k][i[k]] - inter[k][i[k] - 4]);
            }

            // update temporaries
            for (size_t k = 0; k < ELL; ++k)
            {
                tx[k] = inter[k][i[k] - 1];
                ty[k] = inter[k][i[k] - 4];
                sx[k] = 0;
                sy[k] = 0;
            }
        }

        // FINAL M_x MULTIPLICATION
        for (size_t k = 0; k < ELL; ++k)
            for (size_t l = 0; l < ELL; ++l)
                sx[k] = sx[k] + mat[k * ELL + l] * tx[l];

        // FINAL Rho ROTATION
        std::rotate(ty, ty + 1, ty + ELL);

        // FINAL M_y MULTIPLICATION
        for (size_t k = 0; k < ELL; ++k)
            for (size_t l = 0; l < ELL; ++l)
                sy[k] = sy[k] + mat[k * ELL + l] * ty[l];

        // bind output
        constrain(sx[0], 1, out[0]);
    }

    void generate_r1cs_witness()
    {
        std::array<size_t, ELL> i{};
        std::array<Field, ELL> tx{};
        std::array<Field, ELL> ty{};
        std::array<Field, ELL> sx{};
        std::array<Field, ELL> sy{};

        for (size_t j = 0; j < std::min(RATE, ELL); ++j)
            tx[j] = val(in[j][0]);

        for (size_t j = ELL; j < RATE; ++j)
            ty[j - ELL] = val(in[j][0]);

        for (size_t j = 0; j < ROUNDS_N; ++j)
        {
            // ADD CONSTANTS
            for (size_t k = 0; k < ELL; ++k)
            {
                tx[k] = tx[k] + rc[j * ELL + k].first;
                ty[k] = ty[k] + rc[j * ELL + k].second;
            }

            // M_x MULTIPLICATION
            for (size_t k = 0; k < ELL; ++k)
                for (size_t l = 0; l < ELL; ++l)
                    sx[k] = sx[k] + mat[k * ELL + l] * tx[l];

            // Rho ROTATION
            std::rotate(ty.begin(), ty.begin() + 1, ty.end());

            // M_y MULTIPLICATION
            for (size_t k = 0; k < ELL; ++k)
                for (size_t l = 0; l < ELL; ++l)
                    sy[k] = sy[k] + mat[k * ELL + l] * ty[l];

            // PSEUDO-HADAMARD
            for (size_t k = 0; k < ELL; ++k)
            {
                tx[k] = sx[k] + sy[k];
                ty[k] = sy[k] + tx[k];
            }

            // Flystel SBOX
            /*
            Flystel performs the following computations:
                1. x_1 = x_0 - (g(y_0 * y_0) + g_i)
                2. y_1 = y_0 - x_1^(a_i)
                3. x_2 = x_1 + g(y_1 * y_1)
            */
            for (size_t k = 0; k < ELL; ++k)
            {
                /*
                1.  x_1 = x_0 - (g * y_0 * y_0 + g_i)    <==>
                    x_1 - x_0 + g_i = -gy_0 * y_0
                */
                val(inter[k][i[k]]) = tx[k] - (g * ty[k] * ty[k] + g_i);
                ++i[k];

                // 2. y_1 = y_0 - x_1^(a_i)    <==>    (y_0 - y_1)^a = x_1
                // t = x_1^a
                tx[k] = val(inter[k][i[k] - 1]);
                Hash::raise_alpha_inv(tx[k]);
                val(inter[k][i[k]]) = ty[k] - tx[k];
                ++i[k];

                // t^2
                val(inter[k][i[k]]) = tx[k] * tx[k];
                ++i[k];
                // t^4
                val(inter[k][i[k]]) = val(inter[k][i[k] - 1]) * val(inter[k][i[k] - 1]);
                ++i[k];

                // 3. x_2 = x_1 + g(y_1 * y_1)
                val(inter[k][i[k]]) = val(inter[k][i[k] - 4]) +
                                      g * val(inter[k][i[k] - 3]) * val(inter[k][i[k] - 3]);
                ++i[k];
            }
            // update temporaries
            for (size_t k = 0; k < ELL; ++k)
            {
                tx[k] = val(inter[k][i[k] - 1]);
                ty[k] = val(inter[k][i[k] - 4]);
                sx[k] = 0;
                sy[k] = 0;
            }
        }

        // FINAL M_x MULTIPLICATION
        for (size_t k = 0; k < ELL; ++k)
            for (size_t l = 0; l < ELL; ++l)
                sx[k] = sx[k] + mat[k * ELL + l] * tx[l];

        // FINAL Rho ROTATION
        std::rotate(ty.begin(), ty.begin() + 1, ty.end());

        // FINAL M_y MULTIPLICATION
        for (size_t k = 0; k < ELL; ++k)
            for (size_t l = 0; l < ELL; ++l)
                sy[k] = sy[k] + mat[k * ELL + l] * ty[l];

        // bind output
        val(out[0]) = sx[0];
    }
};
