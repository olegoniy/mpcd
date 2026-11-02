#include <algorithm>
#include <utility>

#include <mpcd/VacfAccumulator.hpp>
#include <mpcd/Observables.hpp>

namespace mpcd
{

VacfAccumulator::VacfAccumulator(std::size_t maxLag, std::size_t originStride)
    : maxLag(maxLag),
      originStride(originStride),
      vacfSum(maxLag + 1, 0.0),
      vacfCount(maxLag + 1, 0)
{}


void VacfAccumulator::addOrigin(const System& system, std::size_t step)
{
    VacfOrigin origin;

    origin.startStep = step;
    origin.referenceVelocities = system.v;

    activeOrigins.push_back(std::move(origin));
}


void VacfAccumulator::sample(const System& system, std::size_t currentStep)
{
    for (auto& origin : activeOrigins)
    {
        const std::size_t lag = currentStep - origin.startStep;
        if (lag <= maxLag)
        {
            vacfSum[lag] += solventVelocityAutocorrelation(system, origin.referenceVelocities);
            vacfCount[lag] += 1;
        }
    }

    removeExpiredOrigins(currentStep);
}


bool VacfAccumulator::shouldAddOrigin(std::size_t step) const
{
    return step % originStride == 0;
}


std::vector<double> VacfAccumulator::averagedVacf() const
{
    std::vector<double> result(maxLag + 1, 0.0);

    for (std::size_t lag = 0; lag <= maxLag; ++lag)
        if (vacfCount[lag] > 0)
            result[lag] = vacfSum[lag]/static_cast<double>(vacfCount[lag]);
    return result;
}


const std::vector<std::size_t>& VacfAccumulator::counts() const
{
    return vacfCount;
}


void VacfAccumulator::removeExpiredOrigins(std::size_t currentStep)
{
    activeOrigins.erase(
        std::remove_if(
            activeOrigins.begin(),
            activeOrigins.end(),
            [&](const VacfOrigin& origin)
            {
                return currentStep - origin.startStep > maxLag;
            }
        ),
        activeOrigins.end()
    );
}

}