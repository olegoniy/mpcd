#include <mpcd/MD.hpp>
#include <mpcd/Force.hpp>

namespace mpcd
{

void velocityVerletStep(Polymer& polymer)
{
    Polymer::ParticleMatrix oldForces = polymer.f;
    polymer.r += polymer.dt * polymer.v + (0.5*polymer.dt*polymer.dt/polymer.m)*oldForces;
    // Periodic conditions
    polymer.wrapAllPositionsPBC();
    // Compute forces at new positions
    computeBondForces(polymer);
    polymer.v +=(0.5*polymer.dt/polymer.m) * (oldForces + polymer.f);
}

PolymerSamples runPolymer(Polymer& polymer, std::size_t steps, std::size_t sampleEvery, std::size_t frameEvery, std::size_t bondVectorEvery)
{
    PolymerSamples samples;
    if (sampleEvery > 0)
        samples.reserveObservables(steps/sampleEvery + 1);

    if (frameEvery > 0)
        samples.reserveFrames(steps/frameEvery + 1);
    
    if (bondVectorEvery > 0)
        samples.reserveBondVectors(steps/bondVectorEvery + 1);

    // Initial forces for velocity Verlet
    computeBondForces(polymer);

    // Save initial state
    if (sampleEvery > 0)
        samples.sampleObservables(polymer, 0);

    if (frameEvery > 0)
        samples.sampleFrame(polymer, 0);

    if (bondVectorEvery > 0)
        samples.sampleBondVectors(polymer, 0);
    

    for (std::size_t step = 1; step <= steps; ++step)
    {
        velocityVerletStep(polymer);

        if (sampleEvery > 0 && step%sampleEvery == 0)
            samples.sampleObservables(polymer,step);

        if (frameEvery > 0 && step%frameEvery == 0)
            samples.sampleFrame(polymer, step);

        if (bondVectorEvery > 0 && step % bondVectorEvery == 0)
            samples.sampleBondVectors(polymer, step);
    }

    return samples;
}

}