/**
 * @file IScalarTimeAveragedEquation.cpp
 * @brief Source of scalar time averaged equation interface
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Equations/IScalarTimeAveragedEquation.hpp"

// Project includes
//

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

   MHDVariant IScalarTimeAveragedEquation::updateStoredSolution(const MHDVariant newData, FieldComponents::Spectral::Id, const int i, const int j, const int k)
   {
      // Only update mean on full timestep
      if(this->mTimeFinished)
      {
         T val = incrementTimeAverage(this->mTimeAvg->point(i, j, k), newData, this->time(), this->mTimestep);
         this->mTimeAvg->setPoint(val, i, j, k);
         return val;
      } else
      {
         return noupdateTimeAverage(this->mTimeAvg->point(i,j,k), newData);
      }
   }
}
}
