#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <mpcd/System.hpp>
#include <mpcd/Sampling.hpp>
#include <mpcd/mpcd.hpp>

namespace py = pybind11;

PYBIND11_MODULE(mpcd_cpp, m)
{
    py::class_<mpcd::System>(m, "System")
        .def(
            py::init<
                std::size_t,
                const Eigen::Vector3d&,
                double,
                double,
                double,
                double,
                double,
                unsigned int
            >()
        )

        .def("initPositionsUniform",
             &mpcd::System::initPositionsUniform)

        .def("initVelocitiesMaxwell",
             &mpcd::System::initVelocitiesNormal)

        .def("removeDrift",
             &mpcd::System::removeDrift)

        .def_property_readonly("N",
            [](const mpcd::System& s) {
                return s.N;
            })

        .def_property_readonly("box",
            [](const mpcd::System& s) {
                return s.box;
            })

        .def_property_readonly("r",
            [](mpcd::System& s) {
                py::array array(
                    py::buffer_info(
                        s.r.data(),
                        sizeof(double),
                        py::format_descriptor<double>::format(),
                        2,
                        {s.r.rows(), s.r.cols()},
                        {
                            static_cast<py::ssize_t>(3 * sizeof(double)),
                            static_cast<py::ssize_t>(sizeof(double))
                        }
                    ),
                    py::cast(&s)
                );

                array.attr("setflags")(false);
                return array;
            })

        .def_property_readonly("v",
            [](mpcd::System& s) {
                py::array array(
                    py::buffer_info(
                        s.v.data(),
                        sizeof(double),
                        py::format_descriptor<double>::format(),
                        2,
                        {s.v.rows(), s.v.cols()},
                        {
                            static_cast<py::ssize_t>(3 * sizeof(double)),
                            static_cast<py::ssize_t>(sizeof(double))
                        }
                    ),
                    py::cast(&s)
                );

                array.attr("setflags")(false);
                return array;
            });

    m.def("stream", &mpcd::stream);
    m.def("collide", &mpcd::collide);
    m.def(
        "run_solvent", &mpcd::runSolvent, 
        py::arg("system"), 
        py::arg("steps"), 
        py::arg("sample_period")
    );

    py::class_<mpcd::SolventSamples>(m, "SolventSamples")
    .def_readonly("steps", &mpcd::SolventSamples::steps)
    .def_readonly("temperature", &mpcd::SolventSamples::temperature)
    .def_readonly("kinetic_energy", &mpcd::SolventSamples::kineticEnergy)
    .def_readonly("momentum", &mpcd::SolventSamples::momentum);
}