#pragma once

#include <cstddef>
#include <random>
#include <Eigen/Dense>

namespace mpcd
{

/**
 * @brief Represents a polymer chain used in the MPCD simulation.
 *
 * Stores monomer positions, velocities, forces, simulation box
 * and physical parameters required for molecular dynamics.
 */
class Polymer
{
public:
    using ParticleMatrix = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

    /// Number of monomers in the polymer chain.
    std::size_t N;
    /// Simulation box dimensions.
    Eigen::Vector3d box;
    /// Monomer mass.
    double m;
    /// Thermal energy.
    double kBT;
    /// Harmonic bond spring constant.
    double k;
    /// Equilibrium bond length.
    double r0;
    /// Molecular dynamics time step.
    double dt;
    /// Monomer positions. Each row corresponds to [x, y, z].
    ParticleMatrix r;
    /// Monomer velocities. Each row corresponds to [vx, vy, vz].
    ParticleMatrix v;
    /// Forces acting on the monomers. Each row corresponds to [fx, fy, fz].
    ParticleMatrix f;
    /// Random number generator used for initialization.
    std::mt19937 rng;

    /**
     * @brief Constructs a polymer object and allocates particle data.
     *
     * The constructor initializes all physical parameters and allocates
     * storage for positions, velocities and forces.
     * Positions and velocities are initialized separately.
     *
     * @param nMonomers Number of monomers.
     * @param box Simulation box dimensions.
     * @param dt Molecular dynamics time step.
     * @param bondLength Equilibrium bond length.
     * @param m Monomer mass.
     * @param k Harmonic bond spring constant.
     * @param kBT Thermal energy.
     * @param seed Seed for the random number generator.
     */
    Polymer(
        std::size_t nMonomers,
        const Eigen::Vector3d& box,
        double dt,
        double bondLength,
        double m,
        double k,
        double kBT,
        unsigned int seed
    ):  N(nMonomers),
        box(box),
        m(m),
        kBT(kBT),
        k(k),
        r0(bondLength),
        dt(dt),
        r(nMonomers, 3),
        v(nMonomers, 3),
        f(nMonomers, 3),
        rng(seed)
    {
        f.setZero();
    }

    /**
     * @brief Initializes monomer positions as a random walk.
     */
    void initPositionsRandomWalk();

    /**
     * @brief Initializes monomer positions as a straight chain.
     */
    void initPositionsLinear();

    /**
     * @brief Initializes monomer velocities from a Maxwell-Boltzmann distribution.
     *
     * Each velocity component is sampled from a normal distribution
     * with variance kBT / m.
     */
    void initVelocitiesNormal();

    /**
     * @brief Removes the center-of-mass velocity of the polymer.
     */
    void removeDrift();

    /**
     * @brief Returns the minimum-image displacement vector between two monomers.
     *
     * Periodic boundary conditions are taken into account.
     *
     *  @param idx1 Index of the first monomer.
     * @param idx2 Index of the second monomer.
     */
    Eigen::Vector3d vecDiffPBC(std::size_t idx1, std::size_t idx2) const;

    /**
     * @brief Returns the minimum-image distance between two monomers.
     *
     * Periodic boundary conditions are taken into account.
     *
     * @param idx1 Index of the first monomer.
     * @param idx2 Index of the second monomer.
     * 
     *
     */
    double distancePBC(std::size_t idx1, std::size_t idx2) const;

    /**
     * @brief Wraps Polymer's beads' positions according to PBC 
     */
    void wrapAllPositionsPBC();

private:
    /**
     * @brief Wraps given position according to PBC
     */
    Eigen::Vector3d wrapPositionPBC(const Eigen::Vector3d& position) const;
    /**
     * @brief Generates a random unit vector with isotropic orientation.
     */
    Eigen::Vector3d randomUnitVector();
};

}