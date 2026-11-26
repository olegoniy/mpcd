#include <mpcd/Polymer.hpp>

#include <cmath>
#include <random>

namespace mpcd
{

void Polymer::initPositionsRandomWalk()
{
    r.setZero();

    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::uniform_real_distribution<double> bondScale(0.95, 1.05);

    r(0, 0) = unit(rng) * box(0);
    r(0, 1) = unit(rng) * box(1);
    r(0, 2) = unit(rng) * box(2);

    for (std::size_t i = 1; i < N; ++i)
    {
        Eigen::Vector3d direction = randomUnitVector();
        Eigen::Vector3d step = r0 * direction * bondScale(rng);

        Eigen::Vector3d newPosition = r.row(i - 1).transpose() + step;
        newPosition = wrapPositionPBC(newPosition);
        r.row(i) = newPosition.transpose();
    }
}

void Polymer::initPositionsLinear()
{
    r.setZero();

    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::uniform_real_distribution<double> bondScale(0.95, 1.05);

    r(0, 0) = unit(rng) * box(0);
    r(0, 1) = unit(rng) * box(1);
    r(0, 2) = unit(rng) * box(2);

    Eigen::Vector3d direction = randomUnitVector();

    Eigen::Vector3d step = r0 * direction;

    for (std::size_t i = 1; i < N; ++i)
    {
        Eigen::Vector3d newPosition = r.row(i - 1).transpose() + step * bondScale(rng);;
        newPosition = wrapPositionPBC(newPosition);
        r.row(i) = newPosition.transpose();
    }
}

void Polymer::initVelocitiesNormal()
{
    const double sigma = std::sqrt(kBT / m);
    std::normal_distribution<double> normal(0.0, sigma);

    for (std::size_t i = 0; i < N; ++i)
    {
        v(i, 0) = normal(rng);
        v(i, 1) = normal(rng);
        v(i, 2) = normal(rng);
    }
}

void Polymer::removeDrift()
{
    Eigen::RowVector3d meanVelocity = v.colwise().mean();
    v.rowwise() -= meanVelocity;
}

Eigen::Vector3d Polymer::vecDiffPBC(std::size_t idx1, std::size_t idx2) const
{
    Eigen::Vector3d vecDiff = r.row(idx2).transpose() - r.row(idx1).transpose();
    vecDiff -= (box.array() * (vecDiff.array() / box.array()).round()).matrix();
    return vecDiff;
}

double Polymer::distancePBC(std::size_t idx1, std::size_t idx2) const
{
    return vecDiffPBC(idx1, idx2).norm();
}

void Polymer::wrapAllPositionsPBC()
{
    for (std::size_t i = 0; i < N; ++i)
    {
        Eigen::Vector3d position = r.row(i).transpose();
        position = wrapPositionPBC(position);
        r.row(i) = position.transpose();
    }
}

Eigen::Vector3d Polymer::wrapPositionPBC(const Eigen::Vector3d& position) const
{
    Eigen::Vector3d wrapped = position;

    for (int d = 0; d < 3; ++d)
    {
        wrapped(d) = std::fmod(wrapped(d), box(d));
        if (wrapped(d) < 0.0)
            wrapped(d) += box(d);
    }

    return wrapped;
}

void Polymer::scaleVelocities(double factor)
{
    v *= factor;
}   

void Polymer::setVelocities(const ParticleMatrix& velocities)
{
    if (velocities.rows() != v.rows() || velocities.cols() != 3)
        throw std::invalid_argument("Velocity matrix must have shape (N, 3)"); 

    v = velocities;
}

void Polymer::wrapPositions()
{
    for (Eigen::Index i = 0; i < r.rows(); ++i)
    {
        for (Eigen::Index d = 0; d < 3; ++d)
        {
            const long long crossing = static_cast<long long>(std::floor(r(i, d) / box(d)));

            images(i, d) += crossing;
            r(i, d) -= crossing * box(d);
        }
    }
}

Polymer::ParticleMatrix Polymer::unwrappedPositions() const
{
    ParticleMatrix positions = r;

    for (Eigen::Index i = 0; i < r.rows(); ++i)
    {
        for (Eigen::Index d = 0; d < 3; ++d)
            positions(i, d) += images(i, d) * box(d);
    }

    return positions;
}

Eigen::Vector3d Polymer::randomUnitVector()
{
    std::normal_distribution<double> normal(0.0, 1.0);

    Eigen::Vector3d v(
        normal(rng),
        normal(rng),
        normal(rng)
    );

    return v.normalized();
}

}