/**
 * @file Factory.cpp
 * @brief Implementation of factory for Predictor-Corrector schemes
 */

// System includes
//

// Project includes
//
#include "QuICC/Timestep/Id/ImexPc2.hpp"
#include "Timestep/PredictorCorrector/Factory.hpp"
#include "Timestep/PredictorCorrector/ImExPC2.hpp"
#include "Timestep/PredictorCorrector/Interface.hpp"

namespace QuICC {

namespace Timestep {

namespace PredictorCorrector {

std::shared_ptr<Timestep::Interface> makeInterface(const std::size_t schemeId,
   const MHDFloat time, const Matrix& cfl, const MHDFloat maxError,
   const Timestep::Interface::ScalarEquation_range& scalEq,
   const Timestep::Interface::VectorEquation_range& vectEq,
   Pseudospectral::Coordinator& pseudo)
{
   std::shared_ptr<Timestep::Interface> iface;

   if (schemeId == Id::ImexPc2::id())
   {
      iface = std::make_shared<Interface<ImExPC2>>(time, cfl, maxError, scalEq,
         vectEq, pseudo);
   }

   return iface;
}

} // namespace PredictorCorrector
} // namespace Timestep
} // namespace QuICC
