/**
 * @file IVectorTimeAveragedEquation.cpp
 * @brief Source of vector time averaged equation interface
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Equations/IVectorTimeAveragedEquation.hpp"

// Project includes
//

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

   MHDVariant IVectorTimeAveragedEquation::updateStoredSolution(const MHDVariant newData, FieldComponents::Spectral::Id compId, const int i, const int j, const int k)
   {
      // Only update mean on full timestep
      if(this->mTimeFinished)
      {
         T val = incrementTimeAverage(this->mTimeAvg->comp(compId).point(i, j, k), newData, this->time(), this->mTimestep);
         this->mTimeAvg->rComp(compId).setPoint(val, i, j, k);
         return val;
      } else
      {
         return noupdateTimeAverage(this->mTimeAvg->comp(compId).point(i, j, k), newData);
      }
   }

}
}
