/**
 * @file IVectorTimeAveragedEquation.cpp
 * @brief Source of vector time averaged equation interface
 */

// System includes
//

// Project includes
//
#include "QuICC/Equations/IVectorTimeAveragedEquation.hpp"
#include "QuICC/Equations/AverageSolutionUpdater.hpp"

namespace QuICC {

namespace Equations {

   IVectorTimeAveragedEquation::IVectorTimeAveragedEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IVectorEquation(spEqParams, spScheme, spBackend), mTimeFinished(true), mTimestep(-1.0)
   {
      EquationData::setTime(-4242.0, false);
   }

   IVectorTimeAveragedEquation::~IVectorTimeAveragedEquation()
   {
   }

   void IVectorTimeAveragedEquation::setTime(const MHDFloat time, const bool finished)
   {
      if(this->time() == -4242.0)
      {
         for(auto it = this->spectralRange().first; it != this->spectralRange().second; it++)
         {
            std::visit([&](auto&& p, auto&& t){t->rComp(*it).setData(p->dom(0).perturbation().comp(*it).data());}, this->spUnknown(), this->mTimeAvg);
         }
      }

      this->mTimeFinished = finished;

      if(this->mTimeFinished)
      {
         this->mTimestep = time - this->time();
         EquationData::setTime(time, finished);
      }
   }

   void IVectorTimeAveragedEquation::setUnknown(Framework::Selector::VariantSharedVectorVariable spUnknown)
   {
      IVectorEquation::setUnknown(spUnknown);

      this->mTimeAvg = std::make_shared<Datatypes::VectorField<typename Framework::Selector::ScalarField<T>,FieldComponents::Spectral::Id> >(std::visit([](auto&& p)->auto&&{return p->dom(0).perturbation();},this->spUnknown()));
   }

   void IVectorTimeAveragedEquation::initSolutionUpdater()
   {
      auto range = this->spectralRange();

      for(auto it = range.first; it != range.second; ++it)
      {
         auto spUp = std::make_shared<AverageSolutionUpdater>();
         this->mSolUps.emplace(*it, spUp);
      }
   }

}
}
