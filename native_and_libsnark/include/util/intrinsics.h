#pragma once

#include <inttypes.h>

#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
    #define INTEL_COMPILER
#endif

#if defined(_WIN32) && !defined(__INTEL_LLVM_COMPILER)
    #include <intrin.h>
    #ifndef _bswap64
        #define _bswap64(x) _byteswap_uint64(x)
    #endif
    #ifndef _popcnt16
        #define _popcnt16(x) __popcnt16(x)
    #endif
    #ifndef _popcnt32
        #define _popcnt32(x) __popcnt(x)
    #endif
    #ifndef _popcnt64
        #define _popcnt64(x) __popcnt64(x)
    #endif
#else
    #define __int8 char
    #define __int16 short
    #define __int32 int
    #define __int64 long long
    #define __uint64 unsigned long long

extern "C" inline __uint64 _udiv128(__uint64 hi, __uint64 lo, __uint64 div, __uint64 *rem)
{
    // High bits go in RDX, low bits in RAX, quotient is in RAX, remainder is in RDX
    asm inline( //
        "divq %4"
        : "=d"(hi), "=a"(lo)
        : "d"(hi), "a"(lo), "rm"(div));

    *rem = hi;

    return lo;
}

    #ifndef _WIN32
        #include <x86intrin.h>
        #define _rotr64 _lrotr
        #define _rotl64 _lrotl

extern "C" inline __uint64 _umul128(__uint64 x, __uint64 y, __uint64 *hi)
{
    // x goes in RAX, y goes in RDX, high bits go in RDX, low bits in RAX
    asm inline( //
        "mul %1"
        : "=a"(x), "=d"(y)
        : "a"(x), "d"(y));
    *hi = y;

    return x;
}
    #endif
#endif

#ifdef INTEL_COMPILER
    #undef INTEL_COMPILER
#endif
