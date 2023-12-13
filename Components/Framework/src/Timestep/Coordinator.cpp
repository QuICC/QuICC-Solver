/**
 * @file Coordinator.cpp
 * @brief Implementation of the timestep coordinator
 */

// System includes
//

// Project includes
//
#include "QuICC/Timestep/Coordinator.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "Timestep/RungeKuttaCB/Factory.hpp"
#include "Timestep/PredictorCorrector/Factory.hpp"

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

   void Coordinator::init(const std::size_t schemeId, const MHDFloat time, const Matrix& cfl, const MHDFloat maxError, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, Pseudospectral::Coordinator& pseudo)
   {
      // Try Runge-Kutta CB schemes
      if(!this->mpImpl)
      {
         this->mpImpl = RungeKuttaCB::makeInterface(schemeId, time, cfl, maxError, scalEq, vectEq, pseudo);
      }

      // Try Predictor-Corrector schemes
      if(!this->mpImpl)
      {
         this->mpImpl = PredictorCorrector::makeInterface(schemeId, time, cfl, maxError, scalEq, vectEq, pseudo);
      }

      // No interface was created
      if(!this->mpImpl)
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
