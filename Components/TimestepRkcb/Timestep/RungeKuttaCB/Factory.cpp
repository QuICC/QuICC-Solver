/**
 * @file Factory.cpp
 * @brief Implementation of factory for Runge-Kutta CB schemes
 */

// System includes
//

// Project includes
//
#include "Timestep/RungeKuttaCB/Factory.hpp"
#include "QuICC/Timestep/Id/ImexRkCb2.hpp"
#include "QuICC/Timestep/Id/ImexRkCb3b.hpp"
#include "QuICC/Timestep/Id/ImexRkCb3c.hpp"
#include "QuICC/Timestep/Id/ImexRkCb3d.hpp"
#include "QuICC/Timestep/Id/ImexRkCb3e.hpp"
#include "QuICC/Timestep/Id/ImexRkCb3f.hpp"
#include "QuICC/Timestep/Id/ImexRkCb4.hpp"
#include "Timestep/RungeKuttaCB/ImExRKCB2.hpp"
#include "Timestep/RungeKuttaCB/ImExRKCB3b.hpp"
#include "Timestep/RungeKuttaCB/ImExRKCB3c.hpp"
#include "Timestep/RungeKuttaCB/ImExRKCB3d.hpp"
#include "Timestep/RungeKuttaCB/ImExRKCB3e.hpp"
#include "Timestep/RungeKuttaCB/ImExRKCB3f.hpp"
#include "Timestep/RungeKuttaCB/ImExRKCB4.hpp"
#include "Timestep/RungeKuttaCB/Interface.hpp"

namespace QuICC {

namespace Timestep {

namespace RungeKuttaCB {

std::shared_ptr<Timestep::Interface> makeInterface(const std::size_t schemeId,
   const MHDFloat time, const Matrix& cfl, const MHDFloat maxError,
   const Timestep::Interface::ScalarEquation_range& scalEq,
   const Timestep::Interface::VectorEquation_range& vectEq,
   Pseudospectral::Coordinator& pseudo)
{
   std::shared_ptr<Timestep::Interface> iface;

   if (schemeId == Id::ImexRkCb2::id())
   {
      iface = std::make_shared<Interface<ImExRKCB2>>(time, cfl, maxError,
         scalEq, vectEq, pseudo);
   }
   else if (schemeId == Id::ImexRkCb3b::id())
   {
      iface = std::make_shared<Interface<ImExRKCB3b>>(time, cfl, maxError,
         scalEq, vectEq, pseudo);
   }
   else if (schemeId == Id::ImexRkCb3c::id())
   {
      iface = std::make_shared<Interface<ImExRKCB3c>>(time, cfl, maxError,
         scalEq, vectEq, pseudo);
   }
   else if (schemeId == Id::ImexRkCb3d::id())
   {
      iface = std::make_shared<Interface<ImExRKCB3d>>(time, cfl, maxError,
         scalEq, vectEq, pseudo);
   }
   else if (schemeId == Id::ImexRkCb3e::id())
   {
      iface = std::make_shared<Interface<ImExRKCB3e>>(time, cfl, maxError,
         scalEq, vectEq, pseudo);
   }
   else if (schemeId == Id::ImexRkCb3f::id())
   {
      iface = std::make_shared<Interface<ImExRKCB3f>>(time, cfl, maxError,
         scalEq, vectEq, pseudo);
   }
   else if (schemeId == Id::ImexRkCb4::id())
   {
      iface = std::make_shared<Interface<ImExRKCB4>>(time, cfl, maxError,
         scalEq, vectEq, pseudo);
   }

   return iface;
}

} // namespace RungeKuttaCB
} // namespace Timestep
} // namespace QuICC
