#pragma once

#include "gadget/field_variable.hpp"
#include "gadget/pb_variable_pp.hpp"

template<typename Arion>
class ArionGadget : public GadgetPP<typename Arion::Field>
{
public:
    using Field = typename Arion::Field;
    using Hash = Arion;
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
    using super = GadgetPP<typename Arion::Field>;
    using LC = typename super::LC;

    using super::constrain;
    using super::val;

    static constexpr size_t N = BRANCH_N - 1;

    static constexpr size_t D2_CONSTR = std::__bit_width(Hash::D2); // assuming D2 = 2^k + 1
    static constexpr size_t D2_CONSTR1 = D2_CONSTR - 1;
    static constexpr size_t INTERn_N = D2_CONSTR * ROUNDS_N;
    static constexpr size_t INTERk_N = 5 * ROUNDS_N;

    std::array<std::vector<PbVariablePP<Field>>, BRANCH_N> inter;

public:
    static size_t size() { return INTERn_N + INTERk_N * (BRANCH_N - 1); }

    ArionGadget(libsnark::protoboard<Field> &pb, const BlockVar &in, const DigVar &out,
                const std::string &ap) :
        super{pb, ap}, in{in}, out{out}, inter{}
    {
        for (size_t i = 0; i < INTERn_N; ++i)
            inter[N].emplace_back(pb, FMT(""));

        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < INTERk_N; ++j)
                inter[i].emplace_back(pb, FMT(""));
    }

    void generate_r1cs_constraints()
    {
        std::array<size_t, BRANCH_N> i{};
        std::array<LC, BRANCH_N> t{};
        LC sigma;
        LC g;

        for (size_t i = 0; i < RATE; ++i)
            t[i] = sigma + in[i][0];

        // FIRST CIRCULANT MATRIX
        sigma = t[0];
        g = t[0];
        for (size_t i = 1; i < BRANCH_N; ++i)
            sigma = sigma + t[i];

        t[0] = sigma;
        for (size_t i = 1; i < BRANCH_N; ++i)
            t[0] = t[0] + Hash::circ_mat[i - 1] * t[i];

        for (size_t i = 1; i < BRANCH_N; ++i)
        {
            std::swap(t[i], g);
            t[i] = Hash::circ_mat[BRANCH_N - 1] * t[i];
            t[i] = t[i] + t[i - 1];
            t[i] = t[i] - sigma;
        }

        for (size_t j = 0; j < ROUNDS_N; ++j)
        {
            // GTDS
            // y = x^(1/4097) <==> y^4097 = x
            // y^4096 * y = x
            i[N] += constrain(inter[N][i[N]], inter[N][i[N] + D2_CONSTR1], t[N]);
            // constrain powers of 2
            for (size_t k = 0; k < D2_CONSTR1; ++k)
                i[N] += constrain(inter[N][i[N]], inter[N][i[N]], inter[N][i[N] - 1]);

            // x^d * g(x) + h(x)
            for (size_t k = BRANCH_N - 2; k != (size_t)~0; --k)
            {
                // x^5
                i[k] += constrain(t[k], t[k], inter[k][i[k]]);
                i[k] += constrain(inter[k][i[k] - 1], inter[k][i[k] - 1], inter[k][i[k]]);
                i[k] += constrain(inter[k][i[k] - 1], t[k], inter[k][i[k]]);

                // sigma = sum_{l=k+1}^{BRANCH_N}{x[l] + f[l]}
                sigma = t[k + 1] + inter[k + 1][i[k + 1] - 1];
                for (size_t l = k + 2; l < BRANCH_N; ++l)
                    sigma = sigma + t[l] + inter[l][i[l] - 1];

                i[k] += constrain(sigma, sigma, inter[k][i[k]]);

                // g(x) = s^2 + a1*s + a2
                g = inter[k][i[k] - 1] + Hash::alpha.first * sigma + Hash::alpha.second;
                // h(x) = s^2 + b*s
                sigma = inter[k][i[k] - 1] + Hash::beta1 * sigma;
                // y = x^d * g(x) + h(x) <==> y - h(x) = x^d * g(x)
                i[k] += constrain(inter[k][i[k] - 2], g, inter[k][i[k]] - sigma);
            }

            // CIRCULANT MATRIX
            sigma = inter[0][i[0] - 1];
            for (size_t k = 1; k < BRANCH_N; ++k)
                sigma = sigma + inter[k][i[k] - 1];

            t[0] = sigma;
            for (size_t k = 1; k < BRANCH_N; ++k)
                t[0] = t[0] + Hash::circ_mat[k - 1] * inter[k][i[k] - 1];

            for (size_t k = 1; k < BRANCH_N; ++k)
            {
                t[k] = Hash::circ_mat[BRANCH_N - 1] * inter[k - 1][i[k - 1] - 1];
                t[k] = t[k] + t[k - 1];
                t[k] = t[k] - sigma;
            }

            // ADD CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = t[k] + Hash::round_c[j * BRANCH_N + k];
        }

        constrain(t[0], 1, out[0]);
    }

    void generate_r1cs_witness()
    {
        std::array<size_t, BRANCH_N> i{};
        std::array<Field, BRANCH_N> t{};
        Field sigma;
        Field g;

        for (size_t i = 0; i < RATE; ++i)
            t[i] = val(in[i][0]);

        // FIRST CIRCULANT MATRIX
        sigma = t[0];
        g = t[0];
        for (size_t i = 1; i < BRANCH_N; ++i)
            sigma += t[i];

        t[0] = sigma;
        for (size_t i = 1; i < BRANCH_N; ++i)
            t[0] += Hash::circ_mat[i - 1] * t[i];

        for (size_t i = 1; i < BRANCH_N; ++i)
        {
            std::swap(t[i], g);
            t[i] *= Hash::circ_mat[BRANCH_N - 1];
            t[i] += t[i - 1];
            t[i] -= sigma;
        }

        for (size_t j = 0; j < ROUNDS_N; ++j)
        {
            // GTDS
            // y = x^(1/4097)
            val(inter[N][i[N] + D2_CONSTR1]) = t[N];
            Hash::pow_e(val(inter[N][i[N] + D2_CONSTR1]));
            for (size_t k = 0; k < D2_CONSTR1; ++k)
            {
                val(inter[N][i[N] + (D2_CONSTR1 - k - 1)]) =
                    val(inter[N][i[N] + (D2_CONSTR1 - k)]) * val(inter[N][i[N] + (D2_CONSTR1 - k)]);
            }
            i[N] += D2_CONSTR;

            // x^d * g(x) + h(x)
            for (size_t k = BRANCH_N - 2; k != (size_t)~0; --k)
            {
                // x^5
                val(inter[k][i[k]]) = t[k] * t[k];
                ++i[k];
                val(inter[k][i[k]]) = val(inter[k][i[k] - 1]) * val(inter[k][i[k] - 1]);
                ++i[k];
                val(inter[k][i[k]]) = val(inter[k][i[k] - 1]) * t[k];
                ++i[k];

                // sigma = sum_{l=k+1}^{BRANCH_N}{x[l] + f[l]}
                sigma = t[k + 1] + val(inter[k + 1][i[k + 1] - 1]);
                for (size_t l = k + 2; l < BRANCH_N; ++l)
                    sigma += t[l] + val(inter[l][i[l] - 1]);

                val(inter[k][i[k]]) = sigma * sigma;
                ++i[k];
                // g(x) = s^2 + a1*s + a2
                g = val(inter[k][i[k] - 1]) + Hash::alpha.first * sigma + Hash::alpha.second;
                // h(x) = s^2 + b*s
                sigma = val(inter[k][i[k] - 1]) + Hash::beta1 * sigma;
                // y = x^d * g(x) + h(x)
                val(inter[k][i[k]]) = val(inter[k][i[k] - 2]) * g + sigma;
                ++i[k];
            }
            // CIRCULANT MATRIX
            sigma = val(inter[0][i[0] - 1]);
            for (size_t k = 1; k < BRANCH_N; ++k)
                sigma += val(inter[k][i[k] - 1]);

            t[0] = sigma;
            for (size_t k = 1; k < BRANCH_N; ++k)
                t[0] += Hash::circ_mat[k - 1] * val(inter[k][i[k] - 1]);

            for (size_t k = 1; k < BRANCH_N; ++k)
            {
                t[k] = Hash::circ_mat[BRANCH_N - 1] * val(inter[k - 1][i[k - 1] - 1]);
                t[k] += t[k - 1];
                t[k] -= sigma;
            }

            // ADD CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] += Hash::round_c[j * BRANCH_N + k];
        }

        val(out[0]) = t[0];
    }
};

/*
template<typename FieldT>
class arion_two_to_one_hash_gadget : public ArionGadget<Arion<FieldT, 2, 1>>
{
public:
    using super = ArionGadget<Arion<FieldT, 2, 1>>;
    using Hash = typename super::Hash;
    using DigVar = typename super::DigVar;
    using Field = typename super::Field;
    using BlockVar = typename super::BlockVar;

    arion_two_to_one_hash_gadget(libsnark::protoboard<Field> &pb, const DigVar &x, const DigVar &y,
                                 const DigVar &out, const std::string &ap) :
        super{pb, BlockVar{x, y}, out, ap}
    {}
};
*/