#pragma once

#include <cstddef>
#include <random>
#include <Eigen/Dense>

using ParticleMatrix = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

namespace mpcd
{

/**
 * @brief Represents the MPCD solvent system.
 *
 * Stores solvent particle positions, velocities, simulation parameters,
 * and the random number generator used during initialization and collisions.
 */
class System
{
public:
    /// Number of solvent particles.
    std::size_t N;
    /// Simulation box dimensions.
    Eigen::Vector3d box;
    /// MPCD collision cell size.
    double a;
    /// MPCD collision time step.
    double h;
    /// Solvent particle mass.
    double m;
    /// Thermal energy.
    double kBT;
    /// MPCD collision rotation angle.
    double alpha;
    /// Solvent particle positions. Each row corresponds to [x, y, z].
    ParticleMatrix r;
    /// Solvent particle velocities. Each row corresponds to [vx, vy, vz].
    ParticleMatrix v;
    /// Random number generator used for initialization and MPCD dynamics.
    std::mt19937 rng;

    /**
     * @brief Constructs an MPCD solvent system.
     *
     * @param N Number of solvent particles.
     * @param box Simulation box dimensions.
     * @param a MPCD collision cell size.
     * @param h MPCD collision time step.
     * @param m Solvent particle mass.
     * @param kBT Thermal energy.
     * @param alpha MPCD collision rotation angle.
     * @param seed Seed for the random number generator.
     */
    System(
        std::size_t N,
        const Eigen::Vector3d& box,
        double a,
        double h,
        double m,
        double kBT,
        double alpha,
        unsigned int seed
    ):  N(N),
        box(box),
        a(a),
        h(h),
        m(m),
        kBT(kBT),
        alpha(alpha),
        r(N, 3),
        v(N, 3),
        rng(seed)
    {}

    /**
     * @brief Initializes particle positions uniformly inside the simulation box.
     */
    void initPositionsUniform();

    /**
     * @brief Initializes particle positions uniformly inside a restricted
     * region of the simulation box.
     *
     * Currently used for non-equilibrium initialization where particles
     * occupy only part of the box.
     */
    void initPositionsUniformUpper();

    /**
     * @brief Initializes solvent velocities from a Maxwell-Boltzmann distribution.
     *
     * Each velocity component is sampled from a normal distribution with
     * variance kBT / m.
     */
    void initVelocitiesNormal();

    /**
     * @brief Initializes a non-Maxwellian velocity distribution.
     *
     * Particles are assigned approximately equal speed with random
     * isotropic directions while keeping the kinetic-energy scale
     * consistent with kBT.
     */
    void initVelocitiesNonMaxwell();

    /**
     * @brief Removes the center-of-mass velocity of the solvent.
     *
     * After this operation, the total solvent momentum is approximately zero
     * up to floating-point precision.
     */
    void removeDrift();
};

}