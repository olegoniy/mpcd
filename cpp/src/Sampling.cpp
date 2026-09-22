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

}