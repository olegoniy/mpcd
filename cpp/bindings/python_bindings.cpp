#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include <mpcd/System.hpp>
#include <mpcd/Polymer.hpp> 
#include <mpcd/Sampling.hpp>
#include <mpcd/Observables.hpp>
#include <mpcd/Force.hpp>
#include <mpcd/mpcd.hpp>
#include <mpcd/MD.hpp>

namespace py = pybind11;

PYBIND11_MODULE(mpcd_cpp, m)
{
    // Solvent bindings
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

        .def("initPositionsUniformUpper",
            &mpcd::System::initPositionsUniformUpper)

        .def("initVelocitiesNormal",
             &mpcd::System::initVelocitiesNormal)

        .def("initVelocitiesNonMaxwell",
            &mpcd::System::initVelocitiesNonMaxwell)

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
        
        .def_property_readonly("a",
            [](const mpcd::System& s) {
                return s.a;
            })
        
        .def_property_readonly("h",
            [](const mpcd::System& s) {
                return s.h;
            })

        .def_property_readonly("m",
            [](const mpcd::System& s) {
                return s.m;
            })

        .def_property_readonly("kBT",
            [](const mpcd::System& s) {
                return s.kBT;
            })

        .def_property_readonly("alpha",
            [](const mpcd::System& s) {
                return s.alpha;
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

    py::class_<mpcd::CellList>(m, "CellList")
        .def_readonly("nx", &mpcd::CellList::nx)
        .def_readonly("ny", &mpcd::CellList::ny)
        .def_readonly("nz", &mpcd::CellList::nz)
        .def_readonly("offsets", &mpcd::CellList::offsets)
        .def_readonly("indices", &mpcd::CellList::indices);

    m.def(
        "distributeToCells", &mpcd::distributeToCells,
        py::arg("positions"),
        py::arg("box"),
        py::arg("a")
    );

    py::class_<mpcd::SolventSamples>(m, "SolventSamples")
    .def_readonly("steps", &mpcd::SolventSamples::steps)
    .def_readonly("temperature", &mpcd::SolventSamples::temperature)
    .def_readonly("kineticEnergy", &mpcd::SolventSamples::kineticEnergy)
    .def_readonly("momentum", &mpcd::SolventSamples::momentum)
    .def(py::self += py::self);

    // Polymer bindings
    py::class_<mpcd::Polymer>(m, "Polymer")
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
        >(),
        py::arg("nMonomers"),
        py::arg("box"),
        py::arg("dt"),
        py::arg("bondLength"),
        py::arg("m"),
        py::arg("k"),
        py::arg("kBT"),
        py::arg("seed")
    )

    .def(
        "initPositionsRandomWalk",
        &mpcd::Polymer::initPositionsRandomWalk
    )

    .def(
        "initPositionsLinear",
        &mpcd::Polymer::initPositionsLinear
    )

    .def(
        "initVelocitiesNormal",
        &mpcd::Polymer::initVelocitiesNormal
    )

    .def(
        "removeDrift",
        &mpcd::Polymer::removeDrift
    )

    .def(
        "distancePBC",
        &mpcd::Polymer::distancePBC,
        py::arg("idx1"),
        py::arg("idx2")
    )

    .def(
        "vecDiffPBC",
        &mpcd::Polymer::vecDiffPBC,
        py::arg("idx1"),
        py::arg("idx2")
    )

    .def_readonly(
        "N",
        &mpcd::Polymer::N
    )

    .def_readonly(
        "box",
        &mpcd::Polymer::box
    )

    .def_readonly(
        "m",
        &mpcd::Polymer::m
    )

    .def_readonly(
        "kBT",
        &mpcd::Polymer::kBT
    )

    .def_readonly(
        "k",
        &mpcd::Polymer::k
    )

    .def_readonly(
        "r0",
        &mpcd::Polymer::r0
    )

    .def_readonly(
        "dt",
        &mpcd::Polymer::dt
    )

    .def_property_readonly(
        "r",
        [](mpcd::Polymer& polymer)
        {
            return py::array(
                {
                    static_cast<py::ssize_t>(polymer.N),
                    static_cast<py::ssize_t>(3)
                },
                {
                    static_cast<py::ssize_t>(
                        3 * sizeof(double)
                    ),
                    static_cast<py::ssize_t>(
                        sizeof(double)
                    )
                },
                polymer.r.data(),
                py::cast(&polymer)
            );
        }
    )

    .def_property_readonly(
        "v",
        [](mpcd::Polymer& polymer)
        {
            return py::array(
                {
                    static_cast<py::ssize_t>(polymer.N),
                    static_cast<py::ssize_t>(3)
                },
                {
                    static_cast<py::ssize_t>(
                        3 * sizeof(double)
                    ),
                    static_cast<py::ssize_t>(
                        sizeof(double)
                    )
                },
                polymer.v.data(),
                py::cast(&polymer)
            );
        }
    )

    .def_property_readonly(
        "f",
        [](mpcd::Polymer& polymer)
        {
            return py::array(
                {
                    static_cast<py::ssize_t>(polymer.N),
                    static_cast<py::ssize_t>(3)
                },
                {
                    static_cast<py::ssize_t>(
                        3 * sizeof(double)
                    ),
                    static_cast<py::ssize_t>(
                        sizeof(double)
                    )
                },
                polymer.f.data(),
                py::cast(&polymer)
            );
        }
    );

    py::class_<mpcd::PolymerSamples>(m, "PolymerSamples")
    .def_readonly(
        "steps",
        &mpcd::PolymerSamples::steps
    )
    .def_readonly(
        "kineticEnergy",
        &mpcd::PolymerSamples::kineticEnergy
    )
    .def_readonly(
        "potentialEnergy",
        &mpcd::PolymerSamples::potentialEnergy
    )
    .def_readonly(
        "totalEnergy",
        &mpcd::PolymerSamples::totalEnergy
    )
    .def_readonly(
        "averageBondLength",
        &mpcd::PolymerSamples::averageBondLength
    )
    .def_readonly(
        "endToEndDistance",
        &mpcd::PolymerSamples::endToEndDistance
    )
    .def_readonly(
        "frameSteps",
        &mpcd::PolymerSamples::frameSteps
    )
    .def_readonly(
        "frames",
        &mpcd::PolymerSamples::frames
    )
    .def_readonly(
        "bondVectorSteps",
        &mpcd::PolymerSamples::bondVectorSteps
    )

    .def_readonly(
        "bondVectors",
        &mpcd::PolymerSamples::bondVectors
    )
    .def_readonly(
        "maxBondLength",
        &mpcd::PolymerSamples::maxBondLength
    )

    .def_readonly(
        "momentum",
        &mpcd::PolymerSamples::momentum
    );

    m.def(
        "polymerKineticEnergy",
        &mpcd::polymerKineticEnergy,
        py::arg("polymer")
    );

    m.def(
        "polymerPotentialEnergy",
        &mpcd::polymerPotentialEnergy,
        py::arg("polymer")
    );

    m.def(
        "polymerTotalEnergy",
        &mpcd::polymerTotalEnergy,
        py::arg("polymer")
    );

    m.def(
        "polymerBondLengths",
        &mpcd::polymerBondLengths,
        py::arg("polymer")
    );

    m.def(
        "polymerBondVectors",
        &mpcd::polymerBondVectors,
        py::arg("polymer")
    );

    m.def(
        "runPolymer",
        &mpcd::runPolymer,
        py::arg("polymer"),
        py::arg("steps"),
        py::arg("sample_every"),
        py::arg("frame_every"),
        py::arg("bond_vector_every")
    );
    

    py::class_<mpcd::CoupledSamples>(m, "CoupledSamples")
    .def_readonly(
        "solvent",
        &mpcd::CoupledSamples::solvent
    )
    .def_readonly(
        "polymer",
        &mpcd::CoupledSamples::polymer
    );

    m.def(
    "runCoupled",
    &mpcd::runCoupled,
    py::arg("system"),
    py::arg("polymer"),
    py::arg("steps"),
    py::arg("sample_every"),
    py::arg("frame_every"),
    py::arg("bond_vector_every")
);

}   