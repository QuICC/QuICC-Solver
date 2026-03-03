/**
 * @file Factory.cpp
 * @brief Implementation of factory for Exponential schemes
 */

// System includes
//

// Project includes
//
#include "QuICC/Timestep/Id/ExpEuler.hpp"
#include "Timestep/Exponential/Factory.hpp"
#include "Timestep/Exponential/ExpEuler.hpp"
#include "Timestep/Exponential/Interface.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

std::shared_ptr<Timestep::Interface> makeInterface(const std::size_t schemeId,
   const MHDFloat time, const Matrix& cfl, const MHDFloat maxError,
   const Timestep::Interface::ScalarEquation_range& scalEq,
   const Timestep::Interface::VectorEquation_range& vectEq,
   Pseudospectral::Coordinator& pseudo)
{
   std::shared_ptr<Timestep::Interface> iface;

   if (schemeId == Id::ExpEuler::id())
   {
      iface = std::make_shared<Interface<ExpEuler>>(time, cfl, maxError, scalEq,
         vectEq, pseudo);
   }

   return iface;
}

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
