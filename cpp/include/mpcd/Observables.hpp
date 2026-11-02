#pragma once

#include <mpcd/System.hpp>
#include <mpcd/Polymer.hpp>
#include <mpcd/mpcd.hpp>
#include <Eigen/Dense>

namespace mpcd 
{

double solventKineticEnergy(const System& system, bool removeCOM = false);
double solventTemperature(const System& system, double kB = 1.0, bool removeCOM = true);
Eigen::Vector3d solventMomentum(const System& system);
double coupledTotalEnergy(const System& system, const Polymer& polymer);
double solventVelocityAutocorrelation(const System& system, const ParticleMatrix& referenceVelocities);
double diffusionCoefficientFromVacf(const SolventSamples& samples, double h, std::size_t maxLag);
double theoreticalSolventDiffusion(const System& system);



double polymerKineticEnergy(const Polymer& polymer);
double polymerPotentialEnergy(const Polymer& polymer);
double polymerTotalEnergy(const Polymer& polymer);
std::vector<double> polymerBondLengths(const Polymer& polymer);
std::vector<Eigen::Vector3d> polymerBondVectors(const Polymer& polymer);
double polymerAverageBondLength(const Polymer& polymer);
double polymerEndToEndDistance(const Polymer& polymer);
Eigen::Vector3d polymerMomentum(const Polymer& polymer);
double polymerMaxBondLength(const Polymer& polymer);


std::vector<double> solventRelativeSpeeds(const System& system);

}