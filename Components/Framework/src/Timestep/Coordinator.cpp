/**
 * @file Coordinator.cpp
 * @brief Implementation of the timestep coordinator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Timestep/Coordinator.hpp"

// Project includes
//
#include "QuICC/Tools/Formatter.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Timestep/InterfaceRKCB.hpp"
#include "QuICC/Timestep/Id/ImexRkCb2.hpp"
#include "QuICC/Timestep/ImExRKCB2.hpp"
#include "QuICC/Timestep/Id/ImexRkCb3b.hpp"
#include "QuICC/Timestep/ImExRKCB3b.hpp"
#include "QuICC/Timestep/Id/ImexRkCb3c.hpp"
#include "QuICC/Timestep/ImExRKCB3c.hpp"
#include "QuICC/Timestep/Id/ImexRkCb3d.hpp"
#include "QuICC/Timestep/ImExRKCB3d.hpp"
#include "QuICC/Timestep/Id/ImexRkCb3e.hpp"
#include "QuICC/Timestep/ImExRKCB3e.hpp"
#include "QuICC/Timestep/Id/ImexRkCb3f.hpp"
#include "QuICC/Timestep/ImExRKCB3f.hpp"
#include "QuICC/Timestep/Id/ImexRkCb4.hpp"
#include "QuICC/Timestep/ImExRKCB4.hpp"
#include "QuICC/Timestep/InterfacePC.hpp"
#include "QuICC/Timestep/Id/ImexPc2.hpp"
#include "QuICC/Timestep/ImExPC2.hpp"

namespace QuICC {

namespace Timestep {

   Coordinator::Coordinator()
   {
   }

   MHDFloat Coordinator::time() const
   {
      return mpImpl->time();
   }

   MHDFloat Coordinator::timestep() const
   {
      return mpImpl->timestep();
   }

   void Coordinator::update()
   {
      this->mpImpl->update();
   }

   bool Coordinator::finishedStep() const
   {
      return this->mpImpl->finishedStep();
   }

   void Coordinator::setSolveTime(const std::size_t timeId)
   {
      this->mpImpl->setSolveTime(timeId);
   }

   void Coordinator::getExplicitInput(const std::size_t opId, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
   {
      this->mpImpl->getExplicitInput(opId, scalEq, vectEq, scalVar, vectVar);
   }

   void Coordinator::tuneAdaptive(const MHDFloat time)
   {
      this->mpImpl->tuneAdaptive(time);
   }

   void Coordinator::adaptTimestep(const Matrix& cfl, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      this->mpImpl->adaptTimestep(cfl, scalEq, vectEq);
   }

   void Coordinator::stepForward(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
   {
      this->mpImpl->stepForward(scalEq, vectEq, scalVar, vectVar);
   }

   void Coordinator::init(const std::size_t schemeId, const MHDFloat time, const Matrix& cfl, const MHDFloat maxError, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      if(schemeId == Id::ImexRkCb2::id())
      {
         this->mpImpl = std::make_shared<InterfaceRKCB<ImExRKCB2> >(time, cfl, maxError, scalEq, vectEq);
      }
      else if(schemeId == Id::ImexRkCb3b::id())
      {
         this->mpImpl = std::make_shared<InterfaceRKCB<ImExRKCB3b> >(time, cfl, maxError, scalEq, vectEq);
      }
      else if(schemeId == Id::ImexRkCb3c::id())
      {
         this->mpImpl = std::make_shared<InterfaceRKCB<ImExRKCB3c> >(time, cfl, maxError, scalEq, vectEq);
      }
      else if(schemeId == Id::ImexRkCb3d::id())
      {
         this->mpImpl = std::make_shared<InterfaceRKCB<ImExRKCB3d> >(time, cfl, maxError, scalEq, vectEq);
      }
      else if(schemeId == Id::ImexRkCb3e::id())
      {
         this->mpImpl = std::make_shared<InterfaceRKCB<ImExRKCB3e> >(time, cfl, maxError, scalEq, vectEq);
      }
      else if(schemeId == Id::ImexRkCb3f::id())
      {
         this->mpImpl = std::make_shared<InterfaceRKCB<ImExRKCB3f> >(time, cfl, maxError, scalEq, vectEq);
      }
      else if(schemeId == Id::ImexRkCb4::id())
      {
         this->mpImpl = std::make_shared<InterfaceRKCB<ImExRKCB4> >(time, cfl, maxError, scalEq, vectEq);
      }
      else if(schemeId == Id::ImexPc2::id())
      {
         this->mpImpl = std::make_shared<InterfacePC<ImExPC2> >(time, cfl, maxError, scalEq, vectEq);
      }
      else
      {
         throw std::logic_error("Unknown timestepping scheme requested");
      }
   }

   void Coordinator::printInfo(std::ostream& stream)
   {
      this->mpImpl->printInfo(stream);
   }
}
}
