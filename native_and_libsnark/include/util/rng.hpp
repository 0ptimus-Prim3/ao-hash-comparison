#pragma once

#include "util/intrinsics.h"

#include <array>
#include <random>

static inline uint64_t rand_u64()
{
    static std::mt19937_64 rng{std::random_device{}()};

    return rng();
}
