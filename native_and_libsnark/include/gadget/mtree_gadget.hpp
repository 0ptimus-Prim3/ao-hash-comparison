#pragma once

#include "gadget/digest_variable_pp.hpp"
#include "gadget/field_variable.hpp"
#include "gadget/pb_variable_pp.hpp"
#include "util/array_utils.hpp"

template<size_t height, typename GadHashT>
class MTreeGadget : public GadgetPP<typename GadHashT::Field>
{
public:
    using super = GadgetPP<typename GadHashT::Field>;
    using GadHash = GadHashT;
    using Field = typename GadHash::Field;

    using super::constrain;
    using super::val;

    static constexpr size_t HEIGHT = height;
    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;
    static constexpr size_t DIGEST_SIZE = GadHash::DIGEST_SIZE;
    static constexpr size_t ARITY = GadHash::BLOCK_SIZE / GadHash::DIGEST_SIZE;
    static constexpr bool HASH_ISBOOLEAN = DIGEST_SIZE < DIGEST_VARS;

    using LC = libsnark::linear_combination<Field>;
    using PbVar = PbVariablePP<Field>;
    using DigVar =
        typename std::conditional_t<HASH_ISBOOLEAN, DigestVariablePP<Field>, FieldVariable<Field>>;
    using Protoboard = libsnark::protoboard<Field>;
    using Level = std::array<DigVar, ARITY>;
    using BoolLevel = std::array<PbVar, ARITY - 1>;

private:
    static constexpr size_t HEIGHT1 = HEIGHT - 1;
    static constexpr size_t ARITY1 = ARITY - 1;

    DigVar trans;
    std::vector<Level> other;
    std::vector<DigVar> inter;
    std::vector<Level> children;
    std::vector<GadHash> hash;
    std::vector<BoolLevel> active;

public:
    const DigVar out;

    MTreeGadget(Protoboard &pb, const DigVar &out, const DigVar &trans,
                const std::vector<Level> &other, const std::string &ap) :
        super{pb, ap}, //
        trans{trans},  //
        other{other},  //
        out{out}       //
    {
        for (size_t i = 0; i < HEIGHT1; ++i)
        {
            // inputs for the hash gadget
            children.emplace_back(make_uniform_array<Level>(pb, DIGEST_VARS, FMT("")));
            // active index for next level (i.e. where the previous level was output)
            active.emplace_back(make_uniform_array<BoolLevel>(pb, FMT("")));

            // result of the hash
            inter.emplace_back(pb, DIGEST_VARS, FMT(""));

            // hash gadget
            if (i == HEIGHT1 - 1)
                hash.emplace_back(pb, children[i], out, FMT(""));
            else
                hash.emplace_back(pb, children[i], inter[i], FMT(""));
        }
    }

    void generate_r1cs_constraints()
    {
        // Constraints for layers
        LC sum{0};

        for (size_t j = 0; j < ARITY1; ++j)
        {
            constrain(active[0][j], active[0][j], active[0][j]);
            sum = sum + active[0][j];

            // iterate over all pieces of a single node
            // z = c ? x : y <==> z = xc + y(1 - c) <==> z - y = c(x - y)
            for (size_t k = 0; k < DIGEST_VARS; ++k)
                constrain(active[0][j], trans[k] - other[0][j][k],
                          children[0][j][k] - other[0][j][k]);
        }

        if constexpr (ARITY > 2) // sum = 0|1
            constrain(sum, sum, sum);

        // last branch
        for (size_t k = 0; k < DIGEST_VARS; ++k)
            constrain(1 - sum, trans[k] - other[0][ARITY1][k],
                      children[0][ARITY1][k] - other[0][ARITY1][k]);

        hash[0].generate_r1cs_constraints();

        for (size_t i = 1; i < HEIGHT1; ++i)
        {
            sum = 0;

            for (size_t j = 0; j < ARITY1; ++j)
            {
                constrain(active[i][j], active[i][j], active[i][j]);
                sum = sum + active[i][j];
                // iterate over all pieces of a single node
                // z = c ? x : y <==> z = xc + y(1 - c) <==> z - y = c(x - y)
                for (size_t k = 0; k < DIGEST_VARS; ++k)
                    constrain(active[i][j], inter[i - 1][k] - other[i][j][k],
                              children[i][j][k] - other[i][j][k]);
            }

            if constexpr (ARITY > 2) // sum = 0|1
                constrain(sum, sum, sum);

            // last branch
            for (size_t k = 0; k < DIGEST_VARS; ++k)
                constrain(1 - sum, inter[i - 1][k] - other[i][ARITY1][k],
                          children[i][ARITY1][k] - other[i][ARITY1][k]);

            hash[i].generate_r1cs_constraints();
        }
    }

    void generate_r1cs_witness(size_t idx)
    {
        size_t rem = idx % ARITY;

        for (size_t j = 0; j < ARITY1; ++j)
        {
            val(active[0][j]) = rem == j;

            if (rem == j)
                for (size_t k = 0; k < DIGEST_VARS; ++k)
                    val(children[0][j][k]) = val(trans[k]);
            else
                for (size_t k = 0; k < DIGEST_VARS; ++k)
                    val(children[0][j][k]) = val(other[0][j][k]);
        }

        if (rem == ARITY1)
            for (size_t k = 0; k < DIGEST_VARS; ++k)
                val(children[0][ARITY1][k]) = val(trans[k]);
        else
            for (size_t k = 0; k < DIGEST_VARS; ++k)
                val(children[0][ARITY1][k]) = val(other[0][ARITY1][k]);

        hash[0].generate_r1cs_witness();

        idx /= ARITY;
        for (size_t i = 1; i < HEIGHT1; ++i, idx /= ARITY)
        {
            rem = idx % ARITY;

            for (size_t j = 0; j < ARITY1; ++j)
            {
                val(active[i][j]) = rem == j;

                if (rem == j)
                    for (size_t k = 0; k < DIGEST_VARS; ++k)
                        val(children[i][j][k]) = val(inter[i - 1][k]);
                else
                    for (size_t k = 0; k < DIGEST_VARS; ++k)
                        val(children[i][j][k]) = val(other[i][j][k]);
            }

            if (rem == ARITY1)
                for (size_t k = 0; k < DIGEST_VARS; ++k)
                    val(children[i][ARITY1][k]) = val(inter[i - 1][k]);
            else
                for (size_t k = 0; k < DIGEST_VARS; ++k)
                    val(children[i][ARITY1][k]) = val(other[i][ARITY1][k]);

            hash[i].generate_r1cs_witness();
        }
    }
};
