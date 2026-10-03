#pragma once
#include <mpcd/System.hpp>
#include <mpcd/Sampling.hpp>

namespace mpcd
{
    struct CellList 
    {
        std::size_t nx;
        std::size_t ny;
        std::size_t nz;

        std::vector<std::size_t> offsets;
        std::vector<std::size_t> indices;  
    };
    
    void stream(System& system);
    CellList distributeToCells(const ParticleMatrix& positions, const Eigen::Vector3d& box, double a);
    void collide(System& system);
    SolventSamples runSolvent(System& system, std::size_t steps, std::size_t samplePeriod);

    void collideCoupled(System& system, Polymer& polymer);

    void coupledStep(System& system, Polymer& polymer);

    CoupledSamples runCoupled(System& system, Polymer& polymer, std::size_t steps, std::size_t sampleEvery, std::size_t frameEvery, std::size_t bondVectorEvery);
}
