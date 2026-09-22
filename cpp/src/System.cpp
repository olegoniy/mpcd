#include <mpcd/System.hpp>

namespace mpcd 
{

    void System::initPositionsUniform()
    {
        std::uniform_real_distribution<double> dist(0.0, 0.1); 

        for(std::size_t i = 0; i < N; ++i)
        {
            r(i, 0) = dist(rng) * box(0);
            r(i, 1) = dist(rng) * box(1);
            r(i, 2) = dist(rng) * box(2);
        }
    }

    void System::initVelocitiesNormal()
    {
        std::normal_distribution<double> vel(0.0, sqrt(kBT/m));

        for(std::size_t i = 0; i < N; ++i)
        {
            v(i, 0) = vel(rng);
            v(i, 1) = vel(rng);
            v(i, 2) = vel(rng);
        }
    }

    void System::removeDrift()
    {
        Eigen::RowVector3d v_com = v.colwise().mean();
        v.rowwise() -= v_com;
    }
}   