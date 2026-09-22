#include <mpcd/Observables.hpp>

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

}