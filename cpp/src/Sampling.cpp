#include <mpcd/Sampling.hpp>
#include <mpcd/Observables.hpp>

namespace mpcd {

void SolventSamples::reserve(std::size_t n)
{
    // doesn't change computation, just makes sure memory is reserved in line, not randomly => faster 
    steps.reserve(n);
    temperature.reserve(n);
    kineticEnergy.reserve(n);
    momentum.reserve(n);
}

void SolventSamples::sample(const System& system, std::size_t step)
{
    steps.push_back(step);
    temperature.push_back(solventTemperature(system));
    kineticEnergy.push_back(solventKineticEnergy(system));
    momentum.push_back(solventMomentum(system));
}

SolventSamples& SolventSamples::operator+=(const SolventSamples& other)
{
    std::size_t offset = steps.empty() ? 0 : steps.back();

    for (std::size_t step : other.steps)
    {
        steps.push_back(offset + step);
    }

    temperature.insert(
        temperature.end(),
        other.temperature.begin(),
        other.temperature.end()
    );

    kineticEnergy.insert(
        kineticEnergy.end(),
        other.kineticEnergy.begin(),
        other.kineticEnergy.end()
    );

    momentum.insert(
        momentum.end(),
        other.momentum.begin(),
        other.momentum.end()
    );

    return *this;
}

}