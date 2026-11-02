#pragma once

#include <cstddef>
#include <vector>

#include <mpcd/System.hpp>

namespace mpcd
{

struct VacfOrigin
{
    std::size_t startStep;
    ParticleMatrix referenceVelocities;
};

class VacfAccumulator
{
public:
    VacfAccumulator(std::size_t maxLag, std::size_t originStride);
    void addOrigin(const System& system, std::size_t step);
    void sample(const System& system, std::size_t currentStep);
    bool shouldAddOrigin(std::size_t step) const;

    std::vector<double> averagedVacf() const;

    const std::vector<std::size_t>& counts() const;

private:
    void removeExpiredOrigins(std::size_t currentStep);

    std::size_t maxLag;
    std::size_t originStride;

    std::vector<VacfOrigin> activeOrigins;
    std::vector<double> vacfSum;
    std::vector<std::size_t> vacfCount;         
};

}