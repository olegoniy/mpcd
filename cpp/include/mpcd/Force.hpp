#pragma once

#include <mpcd/Polymer.hpp>

namespace mpcd
{

/**
 * @brief Computes harmonic bond forces acting on the polymer.
 *
 * The force between neighboring monomers is derived from
 *
 * U = 1/2 k (r - r0)^2.
 *
 * Periodic boundary conditions are taken into account.
 *
 * The resulting forces are written to polymer.f.
 */
void computeBondForces(Polymer& polymer);

}