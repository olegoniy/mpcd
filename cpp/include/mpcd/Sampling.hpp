#pragma once

#include <mpcd/System.hpp>
#include <mpcd/Polymer.hpp>

#include <vector>
#include <Eigen/Dense>

namespace mpcd {

struct SolventSamples 
{
    std::vector<std::size_t> steps;
    std::vector<double> temperature;
    std::vector<double> kineticEnergy;
    std::vector<Eigen::Vector3d> momentum;
    std::vector<double> relativeSpeeds;
    std::vector<std::size_t> relativeSpeedSteps;
    std::vector<std::size_t> vacfSteps;
    std::vector<double> velocityAutocorrelation;
    std::vector<std::size_t> averagedVacfLags;
    std::vector<double> averagedVelocityAutocorrelation;
    std::vector<std::size_t> averagedVacfCounts;
    
    

    void reserve(std::size_t n);
    void sample(const System& system, std::size_t step);
    void reserveRelativeSpeeds(std::size_t nSamples, std::size_t nParticles);
    void sampleRelativeSpeeds(const System& system, std::size_t step);
    void reserveVacf(std::size_t n);
    void sampleVacf(const System& system, const ParticleMatrix& referenceVelocities, std::size_t step);



    SolventSamples& operator+=(const SolventSamples& other);
};

struct PolymerSamples
{
    std::vector<std::size_t> steps;

    std::vector<double> kineticEnergy;
    std::vector<double> potentialEnergy;
    std::vector<double> totalEnergy;

    std::vector<double> averageBondLength;
    std::vector<double> endToEndDistance;
    std::vector<std::size_t> bondVectorSteps;
    std::vector<std::vector<Eigen::Vector3d>> bondVectors;

    std::vector<std::size_t> frameSteps;
    std::vector<Polymer::ParticleMatrix> frames;

    void reserve(std::size_t n);

    void reserveObservables(std::size_t n);
    void reserveFrames(std::size_t n);
    void reserveBondVectors(std::size_t n);

    void sampleObservables(const Polymer& polymer, std::size_t step);
    void sampleFrame(const Polymer& polymer, std::size_t step);
    void sampleBondVectors(const Polymer& polymer,std::size_t step);
    std::vector<double> maxBondLength;
    std::vector<Eigen::Vector3d> momentum;
};

struct CoupledSamples
{
    SolventSamples solvent;
    PolymerSamples polymer;
};
}