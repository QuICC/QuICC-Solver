/**
 * @file IFieldEquation.cpp
 * @brief Source of generic field equation interface
 */

// System includes
//

// Project includes
//
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#ifdef QUICC_DEBUG
   #include "QuICC/PhysicalNames/Coordinator.hpp"
   #include "QuICC/SolveTiming/Coordinator.hpp"
   #include "QuICC/Tools/IdToHuman.hpp"
#endif // QUICC_DEBUG

namespace QuICC {

namespace Equations {

   IFieldEquation::IFieldEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IEquation(spEqParams, spScheme, spBackend)
   {
   }

   IFieldEquation::IFieldEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend, std::shared_ptr<EquationOptions> spOptions)
      : IEquation(spEqParams, spScheme, spBackend, spOptions)
   {
   }

   bool IFieldEquation::applyConstraint(FieldComponents::Spectral::Id compId, const std::size_t timeId)
   {
      bool changedSolution = false;

      // Use source kernel
      if(this->mConstraintKernel.count(compId) > 0)
      {
         DebuggerMacro_msg("Apply constraint kernel for " + PhysicalNames::Coordinator::tag(this->name()) + "(" + QuICC::Tools::IdToHuman::toString(compId) + ") at " + SolveTiming::Coordinator::tag(timeId) , 6);

         changedSolution = true;
         this->mConstraintKernel.find(compId)->second->apply(timeId);
      }

      return changedSolution;
   }

   std::vector<std::tuple<MHDVariant,int,int,int>> IFieldEquation::correctionConstraint(FieldComponents::Spectral::Id compId, const std::size_t timeId)
   {
      std::vector<std::tuple<MHDVariant,int,int,int>> corr;
      // Use source kernel
      if(this->mConstraintKernel.count(compId) > 0)
      {
         DebuggerMacro_msg("Correction from constraint kernel for " + PhysicalNames::Coordinator::tag(this->name()) + "(" + QuICC::Tools::IdToHuman::toString(compId) + ") at " + SolveTiming::Coordinator::tag(timeId) , 6);

         corr = this->mConstraintKernel.find(compId)->second->correction(timeId);
      }

      return corr;
   }

   std::shared_ptr<SolutionUpdater> IFieldEquation::spSolutionUpdater(FieldComponents::Spectral::Id compId) const
   {
      if(this->mSolUps.count(compId) == 0)
      {
         throw std::logic_error("Solution updater is not setup properly!");
      }

      return this->mSolUps.at(compId);
   }

   Spectral::Kernel::SharedISpectralKernel IFieldEquation::spSourceKernel(FieldComponents::Spectral::Id compId) const
   {
      // Get source kernel
      if(this->mSrcKernel.count(compId) > 0)
      {
         return this->mSrcKernel.at(compId);
      }
      else
      {
         // This implementation should never get called!
         throw std::logic_error("Requested source kernel without implementation!");
      }
   }

   MHDVariant IFieldEquation::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {
      // Use source kernel
      if(this->mSrcKernel.count(compId) > 0)
      {
         return this->mSrcKernel.find(compId)->second->compute(i, j, k);
      } else
      {
         // This implementation should never get called!
         throw std::logic_error("Activated source term without implementation!");

         return MHDVariant();
      }
   }

   Spectral::Kernel::SharedISpectralKernel IFieldEquation::spBoundaryKernel(FieldComponents::Spectral::Id compId) const
   {
      // Get boundary kernel
      if(this->mBoundaryKernel.count(compId) > 0)
      {
         return this->mBoundaryKernel.at(compId);
      }
      else
      {
         // This implementation should never get called!
         throw std::logic_error("Activated boundary value without implementation!");
      }
   }

   MHDVariant IFieldEquation::boundaryValue(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {
      // Use source kernel
      if(this->mBoundaryKernel.count(compId) > 0)
      {
         return this->mBoundaryKernel.find(compId)->second->compute(i, j, k);
      } else
      {
         // This implementation should never get called!
         throw std::logic_error("Activated boundary value without implementation!");

         return MHDVariant();
      }
   }
} // Equations
} // QuICC
