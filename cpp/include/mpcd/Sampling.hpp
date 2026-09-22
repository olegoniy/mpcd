#pragma once

#include <mpcd/System.hpp>
#include <vector>
#include <Eigen/Dense>

namespace mpcd {

struct SolventSamples {
    std::vector<std::size_t> steps;
    std::vector<double> temperature;
    std::vector<double> kineticEnergy;
    std::vector<Eigen::Vector3d> momentum;

    void reserve(std::size_t n);
    void sample(const System& system, std::size_t step);
};

}