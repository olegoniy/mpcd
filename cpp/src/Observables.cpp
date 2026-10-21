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

double solventVelocityAutocorrelation(const System& system, const ParticleMatrix& referenceVelocities)
{
    if (referenceVelocities.rows() != system.v.rows() || referenceVelocities.cols() != system.v.cols())
        throw std::invalid_argument("Reference velocities have wrong shape");

    double correlation = 0.0;

    for (Eigen::Index i = 0; i < system.v.rows(); ++i)
        correlation += system.v.row(i).dot(referenceVelocities.row(i));

    return correlation/static_cast<double>(system.v.rows());
}

double diffusionCoefficientFromVacf(const SolventSamples& samples, double h, std::size_t maxLag)
{
    if (samples.velocityAutocorrelation.empty())
        throw std::invalid_argument("VACF samples are empty");

    const std::size_t availableLag = samples.velocityAutocorrelation.size() - 1;
    const std::size_t lastLag = std::min(maxLag, availableLag);

    double diffusion = h*samples.velocityAutocorrelation[0]/6.0;

    for (std::size_t lag = 1; lag <= lastLag; ++lag)
        diffusion += h*samples.velocityAutocorrelation[lag]/3.0;

    return diffusion;
}

double theoreticalSolventDiffusion(const System& system)
{
    const double nCells = std::ceil(system.box[0]/system.a)*std::ceil(system.box[1]/system.a)*std::ceil(system.box[2]/system.a);
    const double meanNc = static_cast<double>(system.N)/nCells;
    const double gamma = 2.0*(1.0 - std::cos(system.alpha))/(3.0*meanNc)*(std::exp(-meanNc) + meanNc - 1.0);

    return system.h*system.kBT/system.m*(1.0/gamma - 0.5);
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

Eigen::Vector3d polymerMomentum(const Polymer& polymer)
{
    Eigen::Vector3d momentum = Eigen::Vector3d::Zero();

    for (std::size_t i = 0; i < polymer.N; ++i)
        momentum += polymer.m * polymer.v.row(i).transpose();

    return momentum;
}

double polymerMaxBondLength(const Polymer& polymer)
{
    if (polymer.N < 2)
        return 0.0;

    double maxLength = 0.0;

    for (std::size_t i = 0; i < polymer.N - 1; ++i)
    {
        const double length = polymer.distancePBC(i, i + 1);
        if (length > maxLength)
            maxLength = length;
    }

    return maxLength;
}

double coupledTotalEnergy(const System& system, const Polymer& polymer)
{
    return solventKineticEnergy(system) + polymerTotalEnergy(polymer);
}

std::vector<double> solventRelativeSpeeds(const System& system)
{
    const CellList cells = distributeToCells(system.r, system.box, system.a);

    const std::size_t nCells = cells.nx*cells.ny*cells.nz;

    std::vector<double> relativeSpeeds;
    relativeSpeeds.reserve(static_cast<std::size_t>(system.v.rows()));

    const double thermalVelocity = std::sqrt(system.kBT/system.m);

    for (std::size_t cellId = 0; cellId < nCells; ++cellId)
    {
        const std::size_t begin = cells.offsets[cellId];
        const std::size_t end = cells.offsets[cellId + 1];
        const std::size_t count = end - begin;

        if (count < 2)
            continue;

        Eigen::Vector3d v_com = Eigen::Vector3d::Zero();
        for (std::size_t k = begin; k < end; ++k)
        {
            const std::size_t i = cells.indices[k];
            v_com += system.v.row(i).transpose();
        }

        v_com /= static_cast<double>(count);

        for (std::size_t k = begin; k < end; ++k)
        {
            const std::size_t i = cells.indices[k];
            const Eigen::Vector3d dv = system.v.row(i).transpose() - v_com;
            relativeSpeeds.push_back(dv.norm() / thermalVelocity);
        }
    }

    return relativeSpeeds;
}

}