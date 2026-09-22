#pragma once
#include <mpcd/System.hpp>
#include <mpcd/Sampling.hpp>

namespace mpcd
{
    void stream(System& system);
    void collide(System& system);
    SolventSamples runSolvent(System& system, std::size_t steps, std::size_t samplePeriod);
}