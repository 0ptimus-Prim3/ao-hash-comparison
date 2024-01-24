#pragma once

#include "intrinsics.h"
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>

template<typename T>
void consume(T &x)
{
    void *p = (T *)(&x);
    // do nothing
    asm volatile("" : : "rm"(p));
}

template<typename Fn> //
double measure(Fn foo, size_t repeat = 1, size_t times = 1, const char *name = nullptr,
               bool output = true)
{
    using clk = std::chrono::high_resolution_clock;

    uint64_t clocks_avg = 0;
    double elap_avg = 0, cpi_avg = 0, ops_avg = 0;
    std::ostream out{std::cout.rdbuf()};

    out << std::left;
    if (output)
        out << "==== Measuring " << repeat << " iteration" << (repeat > 1 ? "s" : "")
            << (name ? " of " : "") << (name ? name : "") << " ====\n"
            << "Clocks          Millisec        CPI             IPS\n";

    for (size_t i = 0; i < times; ++i)
    {
        uint64_t c_start = __rdtsc();
        auto start = clk::now();
        for (size_t j = 0; j < repeat; ++j)
        {
            if constexpr (std::is_invocable_v<Fn, size_t>)
                foo(j);
            else
                foo();
        }

        uint64_t clocks = __rdtsc() - c_start;
        uint64_t elap = (clk::now() - start).count();

        double elap_ms = elap / 1'000'000.;
        double cpi = (double)clocks / repeat;
        double ops = repeat * 1'000'000'000. / elap;

        clocks_avg = (clocks_avg * i + clocks) / (i + 1);
        elap_avg = (elap_avg * i + elap) / (i + 1);
        cpi_avg = (cpi_avg * i + cpi) / (i + 1);
        ops_avg = (ops_avg * i + ops) / (i + 1);
        if (output)
            out << std::setw(16) << clocks << std::setw(16) << elap_ms << std::setw(16) << cpi
                << ops << '\n';
    }

    elap_avg /= 1'000'000;
    if (output)
    {
        if (times > 1)
            out << "---- Average ----\n"
                << std::setw(16) << clocks_avg << std::setw(16) << elap_avg << std::setw(16)
                << cpi_avg << std::setw(16) << ops_avg << '\n';
        out << "\n\n";
    }

    return elap_avg;
}
