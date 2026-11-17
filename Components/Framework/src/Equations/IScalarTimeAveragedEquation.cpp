/**
 * @file IScalarTimeAveragedEquation.cpp
 * @brief Source of scalar time averaged equation interface
 */

// System includes
//

// Project includes
//
#include "QuICC/Equations/IScalarTimeAveragedEquation.hpp"
#include "QuICC/Equations/AverageSolutionUpdater.hpp"

namespace QuICC {

namespace Equations {

   IScalarTimeAveragedEquation::IScalarTimeAveragedEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IScalarEquation(spEqParams, spScheme, spBackend), mTimeFinished(false), mTimestep(-1.0)
   {
      EquationData::setTime(-4242.0, false);
   }

   IScalarTimeAveragedEquation::~IScalarTimeAveragedEquation()
   {
   }

   void IScalarTimeAveragedEquation::setTime(const MHDFloat time, const bool finished)
   {
      if(this->time() == -4242.0)
      {
         std::visit([&](auto&& p, auto&& t){t->setData(p->dom(0).perturbation().data());}, this->spUnknown(), this->mTimeAvg);
         EquationData::setTime(time, finished);
      }

      this->mTimeFinished = finished;

      if(this->mTimeFinished)
      {
         this->mTimestep = time - this->time();
         EquationData::setTime(time, finished);
      }
   }

   void IScalarTimeAveragedEquation::setUnknown(Framework::Selector::VariantSharedScalarVariable spUnknown)
   {
      IScalarEquation::setUnknown(spUnknown);

      this->mTimeAvg = std::make_shared<typename Framework::Selector::ScalarField<T> >(std::visit([](auto&& p)->auto&&{return p->dom(0).perturbation();},this->spUnknown()));
   }

   void IScalarTimeAveragedEquation::initSolutionUpdater()
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
