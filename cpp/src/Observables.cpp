#include <mpcd/Observables.hpp>
#include <vector>

namespace mpcd 
{
double solventKineticEnergy(const System& system, bool removeCOM)
{
    Eigen::RowVector3d vCOM{0,0,0};
    if (removeCOM) 
        vCOM = system.v.colwise().mean();

    return 0.5 * system.m * (system.v.rowwise() - vCOM).rowwise().squaredNorm().sum();
}

double solventTemperature(const System& system, double kB, bool removeCOM)
{
    std::size_t dof = 3 * system.N;
    if (removeCOM)
        dof -= 3;
    double K =solventKineticEnergy(system, removeCOM);

    return 2*K/(dof * kB);
}

Eigen::Vector3d solventMomentum(const System& system)
{
    return system.m * system.v.colwise().sum().transpose();
}

double polymerKineticEnergy(const Polymer &polymer)
{
    return 0.5 * polymer.m * polymer.v.array().square().sum();
}

double polymerPotentialEnergy(const Polymer& polymer)
{
    if (polymer.N < 2)
        return 0.0;

    double potentialEnergy = 0.0;
    for (std::size_t i = 0; i < polymer.N - 1; ++i)
    {
        double bondLength = polymer.distancePBC(i, i + 1);
        double extension = bondLength - polymer.r0;

        potentialEnergy += 0.5 * polymer.k * extension * extension;
    }

    return potentialEnergy;
}

double polymerTotalEnergy(const Polymer& polymer)
{
    return polymerKineticEnergy(polymer) + polymerPotentialEnergy(polymer);
}

std::vector<double> polymerBondLengths(const Polymer& polymer)
{
    std::vector<double> lengths;
    if (polymer.N < 2)
        return lengths;

    lengths.reserve(polymer.N - 1);

    for (std::size_t i = 0; i < polymer.N - 1; ++i)
        lengths.push_back(polymer.distancePBC(i, i + 1));

    return lengths;
}

std::vector<Eigen::Vector3d> polymerBondVectors(const Polymer& polymer)
{
    std::vector<Eigen::Vector3d> vectors;
    if (polymer.N < 2)
        return vectors;

    vectors.reserve(polymer.N - 1);

    for (std::size_t i = 0; i < polymer.N - 1; ++i)
        vectors.push_back(polymer.vecDiffPBC(i, i + 1));

    return vectors;
}

double polymerAverageBondLength(const Polymer& polymer)
{
    if (polymer.N < 2)
        return 0.0;

    double sum = 0.0;
    for (std::size_t i = 0; i < polymer.N - 1; ++i)
        sum += polymer.distancePBC(i, i + 1);

    return sum/static_cast<double>(polymer.N - 1);
}

double polymerEndToEndDistance(const Polymer& polymer)
{
    if (polymer.N < 2)
        return 0.0;

    return polymer.distancePBC(0, polymer.N - 1);
}

}