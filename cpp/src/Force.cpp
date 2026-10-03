#include <mpcd/Force.hpp>

namespace mpcd
{

float bondForce(double k, double distance, double bondLength)
{
    return k*(distance - bondLength);
}

void computeBondForces(Polymer& polymer)
{
    polymer.f.setZero();

    for (std::size_t i = 0; i < polymer.N - 1; ++i)
    {
        Eigen::Vector3d bond = polymer.vecDiffPBC(i, i + 1);
        double distance = bond.norm();
        Eigen::Vector3d direction = bond / distance;
        Eigen::Vector3d force = bondForce(polymer.k, distance, polymer.r0)*direction;

        polymer.f.row(i) += force.transpose();
        polymer.f.row(i + 1) -= force.transpose();
    }
}

}