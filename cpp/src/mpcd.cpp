#include <mpcd/mpcd.hpp>

#include <mpcd/System.hpp>
#include <mpcd/MD.hpp>
#include <mpcd/Polymer.hpp>
#include <mpcd/Sampling.hpp>
#include <mpcd/Force.hpp>

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace {

    void wrapPositions(ParticleMatrix& positions, const Eigen::Vector3d& box)
    {
        for (Eigen::Index i = 0; i < positions.rows(); ++i) 
        {
            for (Eigen::Index d = 0; d < 3; ++d) 
            {
                positions(i, d) = std::fmod(positions(i, d), box(d));

                if (positions(i, d) < 0.0)
                    positions(i, d) += box(d);
            }
        }
    }

    Eigen::Matrix3d generateRotation(std::mt19937& rng, double alpha)
    {
        std::uniform_real_distribution<double> phiDist(0.0, 2.0 * std::numbers::pi);

        std::uniform_real_distribution<double> zDist(-1.0, 1.0);

        const double phi = phiDist(rng);
        const double z = zDist(rng);

        const double r = std::sqrt(1.0 - z * z);

        const double Rx = r * std::cos(phi);
        const double Ry = r * std::sin(phi);
        const double Rz = z;

        const double c = std::cos(alpha);
        const double s = std::sin(alpha);

        Eigen::Matrix3d R;

        R <<
            Rx*Rx + (1.0 - Rx*Rx)*c,
            Rx*Ry*(1 - c) - Rz*s,
            Rx*Rz*(1 - c) + Ry*s,

            Rx*Ry*(1 - c) + Rz*s,
            Ry*Ry + (1.0 - Ry*Ry)*c,
            Ry*Rz*(1 - c) - Rx*s,

            Rx*Rz*(1 - c) - Ry*s,
            Ry*Rz*(1 - c) + Rx*s,
            Rz*Rz + (1.0 - Rz*Rz)*c;

        return R;
    }

    void rotateInCell(mpcd::System& system, const mpcd::CellList& cells, std::size_t cellId, const Eigen::Matrix3d& rotation)
    {
        const std::size_t begin = cells.offsets[cellId];
        const std::size_t end = cells.offsets[cellId + 1];

        const std::size_t count = end - begin;

        if (count < 2) // no need to rotate, if there is 1 or 0 particles
            return;

        Eigen::Vector3d v_com = Eigen::Vector3d::Zero();

        for (std::size_t k = begin; k < end; ++k) {
            const std::size_t i = cells.indices[k];
            v_com += system.v.row(i).transpose();
        }

        v_com /= static_cast<double>(count);

        for (std::size_t k = begin; k < end; ++k) {
            const std::size_t i = cells.indices[k];

            Eigen::Vector3d dv =
                system.v.row(i).transpose() - v_com;

            system.v.row(i) =
                (v_com + rotation * dv).transpose();
        }
    }

    Eigen::Vector3d coupledCMVelocity(const mpcd::System& system, const mpcd::CellList& solventCells, const mpcd::Polymer& polymer, const mpcd::CellList& polymerCells, std::size_t cellId)
    {
        const std::size_t solventBegin = solventCells.offsets[cellId];
        const std::size_t solventEnd = solventCells.offsets[cellId + 1];

        const std::size_t polymerBegin = polymerCells.offsets[cellId];
        const std::size_t polymerEnd = polymerCells.offsets[cellId + 1];

        Eigen::Vector3d momentum = Eigen::Vector3d::Zero();
        double totalMass = 0.0;

        for (std::size_t k = solventBegin; k < solventEnd; ++k) 
        {
            const std::size_t i = solventCells.indices[k];
            momentum += system.m*system.v.row(i).transpose();
            totalMass += system.m;
        }

        for (std::size_t k = polymerBegin; k < polymerEnd; ++k)
        {
            const std::size_t i = polymerCells.indices[k];
            momentum += polymer.m*polymer.v.row(i).transpose();
            totalMass += polymer.m;
        }

        return momentum / totalMass;
    }

    void rotateCoupledCell(mpcd::System& system, const mpcd::CellList& solventCells, mpcd::Polymer& polymer, const mpcd::CellList& polymerCells, std::size_t cellId, const Eigen::Matrix3d& rotation)
    {
        const std::size_t solventBegin = solventCells.offsets[cellId];
        const std::size_t solventEnd = solventCells.offsets[cellId + 1];

        const std::size_t polymerBegin = polymerCells.offsets[cellId];
        const std::size_t polymerEnd = polymerCells.offsets[cellId + 1];

        const Eigen::Vector3d v_com = coupledCMVelocity(system, solventCells, polymer, polymerCells, cellId);

        for (std::size_t k = solventBegin; k < solventEnd; ++k)
        {
            const std::size_t i = solventCells.indices[k];
            Eigen::Vector3d dv = system.v.row(i).transpose() - v_com;
            system.v.row(i) = (v_com + rotation * dv).transpose();
        }

        for (std::size_t k = polymerBegin; k < polymerEnd; ++k)
        {
            const std::size_t i = polymerCells.indices[k];
            Eigen::Vector3d dv = polymer.v.row(i).transpose() - v_com;
            polymer.v.row(i) = (v_com + rotation * dv).transpose();
        }
    }
}

namespace mpcd 
{
    void stream(System& system)
    {
        system.r += system.h * system.v;
        wrapPositions(system.r, system.box);
    }

    CellList distributeToCells(const ParticleMatrix& positions, const Eigen::Vector3d& box, double a)
    {
        // Create correct sized grid
        const std::size_t nx = static_cast<std::size_t>(std::ceil(box(0) / a));
        const std::size_t ny = static_cast<std::size_t>(std::ceil(box(1) / a));
        const std::size_t nz = static_cast<std::size_t>(std::ceil(box(2) / a));

        const std::size_t nCells = nx * ny * nz;
        const std::size_t N = static_cast<std::size_t>(positions.rows());

        CellList cells{
            nx,
            ny,
            nz,
            std::vector<std::size_t>(nCells + 1, 0),
            std::vector<std::size_t>(N)
        };

        // actual distribution to cells (numbers at first)
        for (Eigen::Index i = 0; i < positions.rows(); ++i)
        {
            const std::size_t ix = static_cast<std::size_t>(positions(i, 0) / a);
            const std::size_t iy = static_cast<std::size_t>(positions(i, 1) / a);
            const std::size_t iz = static_cast<std::size_t>(positions(i, 2) / a);

            const std::size_t cellId = ix + nx * (iy + ny * iz);
            ++cells.offsets[cellId + 1];
        }

        // convert numbers into offsets
        for (std::size_t c = 0; c < nCells; ++c) {
            cells.offsets[c + 1] += cells.offsets[c];
        }

        // create oneline particles' indicies list
        std::vector<std::size_t> cursor = cells.offsets;

        for (Eigen::Index i = 0; i < positions.rows(); ++i) {

            const std::size_t ix = static_cast<std::size_t>(positions(i, 0) / a);
            const std::size_t iy = static_cast<std::size_t>(positions(i, 1) / a);
            const std::size_t iz = static_cast<std::size_t>(positions(i, 2) / a);

            const std::size_t cellId = ix + nx * (iy + ny * iz);

            cells.indices[cursor[cellId]] = static_cast<std::size_t>(i);
            ++cursor[cellId];
        }

        return cells;
    }

    void collide(System& system)
    {
        // random grid shift 
        std::uniform_real_distribution<double> dist(-system.a/2, system.a/2);
        Eigen::Vector3d shift(
            dist(system.rng),
            dist(system.rng),
            dist(system.rng)
        );

        ParticleMatrix shifted_positions = system.r.rowwise() + shift.transpose();        
        // periodic wraping  
        wrapPositions(shifted_positions, system.box);
        //distribution to cells
        CellList cells = distributeToCells(shifted_positions, system.box, system.a);

        const std::size_t nCells = cells.nx * cells.ny * cells.nz;
        for (std::size_t cellId = 0; cellId < nCells; ++cellId) 
        {

            const std::size_t begin = cells.offsets[cellId];
            const std::size_t end   = cells.offsets[cellId + 1];
            const std::size_t count = end - begin;

            if (count > 1) 
            {
                const Eigen::Matrix3d rotation = generateRotation(system.rng, system.alpha);
                rotateInCell(system, cells, cellId, rotation);
            }
        }
    }

    void collideCoupled(System& system, Polymer& polymer)
    {
        std::uniform_real_distribution<double> dist(-system.a/2, system.a/2);

        Eigen::Vector3d shift(
            dist(system.rng),
            dist(system.rng),
            dist(system.rng)
        );

        ParticleMatrix shiftedSolvent = system.r.rowwise() + shift.transpose();
        Polymer::ParticleMatrix shiftedPolymer = polymer.r.rowwise() + shift.transpose();

        wrapPositions(shiftedSolvent, system.box);
        wrapPositions(shiftedPolymer, polymer.box);

        CellList solventCells = distributeToCells(shiftedSolvent, system.box, system.a);
        CellList polymerCells = distributeToCells(shiftedPolymer, polymer.box, system.a);

        const std::size_t nCells = solventCells.nx * solventCells.ny * solventCells.nz;

        for (std::size_t cellId = 0; cellId < nCells; ++cellId)
        {
            const std::size_t solventBegin = solventCells.offsets[cellId];
            const std::size_t solventEnd = solventCells.offsets[cellId + 1];

            const std::size_t polymerBegin = polymerCells.offsets[cellId];
            const std::size_t polymerEnd = polymerCells.offsets[cellId + 1];

            const std::size_t solventCount = solventEnd - solventBegin;
            const std::size_t polymerCount = polymerEnd - polymerBegin;

            if (solventCount + polymerCount > 1)
            {
                const Eigen::Matrix3d rotation = generateRotation(system.rng, system.alpha);
                rotateCoupledCell(system, solventCells, polymer, polymerCells,cellId, rotation);
            }
        }
    }

    SolventSamples runSolvent(System& system, std::size_t steps, std::size_t samplePeriod)
    {
        SolventSamples samples;

        if (samplePeriod > 0) {
            samples.reserve(steps / samplePeriod);
        }

        for (std::size_t step = 0; step < steps; ++step)
        {
            stream(system);
            collide(system);
            if (samplePeriod > 0 && (step + 1) % samplePeriod == 0)
                samples.sample(system, step + 1);
        }

        return samples;
    }

    void coupledStep(System& system, Polymer& polymer)
    {
        const double ratio = system.h / polymer.dt;
        const std::size_t mdSteps = static_cast<std::size_t>(std::llround(ratio));

        if (std::abs(ratio - static_cast<double>(mdSteps)) > 1e-10)
            throw std::runtime_error("system.h must be an integer multiple of polymer.dt");

        for (std::size_t step = 0; step < mdSteps; ++step)
            velocityVerletStep(polymer);

        stream(system);
        collideCoupled(system, polymer);
    }

    CoupledSamples runCoupled(System& system, Polymer& polymer, std::size_t steps, std::size_t sampleEvery, std::size_t frameEvery, std::size_t bondVectorEvery)
    {
        CoupledSamples samples;

        if (sampleEvery > 0)
        {
            samples.solvent.reserve(steps / sampleEvery + 1);
            samples.polymer.reserveObservables(steps / sampleEvery + 1);
        }

        if (frameEvery > 0)
            samples.polymer.reserveFrames(steps / frameEvery + 1);

        if (bondVectorEvery > 0)
            samples.polymer.reserveBondVectors(steps / bondVectorEvery + 1);

        computeBondForces(polymer);

        if (sampleEvery > 0)
        {
            samples.solvent.sample(system, 0);
            samples.polymer.sampleObservables(polymer, 0);
        }

        if (frameEvery > 0)
            samples.polymer.sampleFrame(polymer, 0);

        if (bondVectorEvery > 0)
            samples.polymer.sampleBondVectors(polymer, 0);

        for (std::size_t step = 1; step <= steps; ++step)
        {
            coupledStep(system, polymer);

            if (sampleEvery > 0 && step % sampleEvery == 0)
            {
                samples.solvent.sample(system, step);
                samples.polymer.sampleObservables(polymer, step);
            }

            if (frameEvery > 0 && step % frameEvery == 0)
                samples.polymer.sampleFrame(polymer, step);

            if (bondVectorEvery > 0 && step % bondVectorEvery == 0)
                samples.polymer.sampleBondVectors(polymer, step);
        }

        return samples;
    }
}