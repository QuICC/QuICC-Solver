/**
 * @file Factory.cpp
 * @brief Implementation of factory for Predictor-Corrector schemes
 */

// System includes
//

// Project includes
//
#include "QuICC/Timestep/Id/ImexPc2.hpp"
#include "Timestep/PredictorCorrector/FactoryViews.hpp"
#include "Timestep/PredictorCorrector/ImExPC2.hpp"
#include "Timestep/PredictorCorrector/InterfaceViews.hpp"

namespace QuICC {

namespace Timestep {

namespace PredictorCorrector {

std::shared_ptr<Timestep::Interface> makeInterfaceViews(const std::size_t schemeId,
   const MHDFloat time, const Matrix& cfl, const MHDFloat maxError,
   const Timestep::Interface::ScalarEquation_range& scalEq,
   const Timestep::Interface::VectorEquation_range& vectEq,
   Pseudospectral::Coordinator& pseudo)
{
   std::shared_ptr<Timestep::Interface> iface;

   if (schemeId == Id::ImexPc2::id())
   {
      iface = std::make_shared<InterfaceViews<ImExPC2>>(time, cfl, maxError, scalEq,
         vectEq, pseudo);
   }

   return iface;
}

} // namespace PredictorCorrector
} // namespace TimestepViews
} // namespace QuICC
