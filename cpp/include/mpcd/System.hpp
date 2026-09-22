#pragma once
#include <cstddef>
#include <vector>
#include <random>
#include <Eigen/Dense> 

using ParticleMatrix = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

namespace mpcd {
    class System {
    public:

    // parameters
        size_t N;
        Eigen::Vector3d box;
        double a;
        double h;
        double m;
        double kBT;
        double alpha;
    // particles data

    ParticleMatrix r;
    ParticleMatrix v;
    // RNG
        std::mt19937 rng;
    // constructor
        System(
            std::size_t N,
            const Eigen::Vector3d& box,
            float a,
            float h,
            float m,
            double kBT,
            double alpha,
            unsigned int seed
        )
        : N(N),
        box(box),
        a(a),
        h(h),
        m(m),
        kBT(kBT),
        alpha(alpha),
        r(N, 3),
        v(N, 3),
        rng(seed)
        {};

        void initPositionsUniform();
        void initVelocitiesNormal();
        void removeDrift();
    };
}