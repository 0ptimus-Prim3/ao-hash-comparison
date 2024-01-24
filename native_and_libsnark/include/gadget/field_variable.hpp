#pragma once

#include "gadget/gadget_pp.hpp"
#include "util/algebra.hpp"
#include "util/string_utils.hpp"
#include <gmpxx.h>

template<typename FieldT>
class FieldVariable : public GadgetPP<FieldT>
{
private:
    static constexpr size_t FIELD_SIZE = field_size<FieldT>();

    libsnark::pb_variable_array<FieldT> vars;

public:
    using super = GadgetPP<FieldT>;

    using super::constrain;
    using super::val;

    FieldVariable(libsnark::protoboard<FieldT> &pb, const size_t num_vars,
                  const std::string &annotation_prefix) :
        super{pb, annotation_prefix}
    {
        vars.allocate(pb, num_vars, annotation_prefix);
    }

    void generate_r1cs_constraints(){};


    template<typename Range>
    void generate_r1cs_witness(const Range &contents)
    {
        auto it = std::begin(contents);
        auto end = std::end(contents);

        if constexpr (std::is_same_v<FieldT, std::remove_reference_t<decltype(*it)>>)
        {
            for (size_t i = 0; i < vars.size() && it != end; ++i, ++it)
                val(vars[i]) = *it;
        }
        else
        {
            static constexpr size_t IT_SIZE = sizeof(decltype(*it));
            uint8_t buff[FIELD_SIZE];

            for (size_t i = 0; i < vars.size() && it != end; ++i)
            {
                for (size_t j = 0; j < FIELD_SIZE; j += IT_SIZE, ++it)
                    memcpy(buff + j, &*it, IT_SIZE);

                val(vars[i]) = field_load<FieldT>(buff);
            }
        }
    }

    void generate_r1cs_witness(const void *data, size_t sz)
    {
        const char *d = (const char *)data;

        for (size_t i = 0, n = std::min(sz / FIELD_SIZE, vars.size()); i < n; ++i)
            val(vars[i]) = field_load<FieldT>(d + i * FIELD_SIZE);
    }

    const auto &operator[](size_t i) const { return vars[i]; }

    auto &operator[](size_t i) { return vars[i]; }

    size_t size() const { return vars.size(); }

    auto begin() { return vars.begin(); }

    auto end() { return vars.end(); }

    auto begin() const { return vars.begin(); }

    auto end() const { return vars.end(); }
};
