#pragma once

#include <mpcd/Polymer.hpp>
#include <mpcd/Sampling.hpp>

namespace mpcd
{

/**
 * @brief Advances the polymer by one molecular-dynamics step
 * using the velocity Verlet algorithm.
 *
 * Bond forces are recomputed after updating the positions.
 */
void velocityVerletStep(Polymer& polymer);
/**
 * @brief Run simulation on a free polymer with sampling
 * 
 * @param steps number of md steps
 * @param sampleEvery number of md steps between sampling
 * @param frameEvery number of md steps to sample all coordinates
 */
PolymerSamples runPolymer(Polymer& polymer, std::size_t steps, std::size_t sampleEvery, std::size_t frameEvery, std::size_t bondVectorEvery);

}