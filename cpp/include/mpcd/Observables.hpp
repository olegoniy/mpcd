#pragma once

#include <mpcd/System.hpp>
#include <Eigen/Dense>

namespace mpcd 
{

double solventKineticEnergy(const System& system, bool removeCOM = false);
double solventTemperature(const System& system, double kB = 1.0, bool removeCOM = true);
Eigen::Vector3d solventMomentum(const System& system);

}