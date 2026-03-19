/**
 * @file ExpEuler.cpp
 * @brief Implementation of an exponential Euler scheme of
 * order 1
 */

// System includes
//

// Project includes
//
#include "Timestep/Exponential/ExpEuler.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

ExpEuler::ExpEuler()
{
}

// Scheme requires 1 substep
int ExpEuler::steps() const
{
   return 1;
}

// Scheme order
int ExpEuler::order() const
{
   return 2;
}

// Scheme has embedded lower order scheme?
bool ExpEuler::hasEmbedded() const
{
   return false;
}

// Name of the scheme
std::string ExpEuler::name() const
{
   return "ExpEuler";
}

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
