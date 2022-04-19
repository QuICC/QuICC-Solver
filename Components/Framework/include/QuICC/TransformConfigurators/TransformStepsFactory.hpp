/**
 * @file TransformStepsFactory.hpp
 * @brief Factory to dispatch to implementation of transform steps depending on spatial scheme
 */

#ifndef QUICC_TRANSFORM_TRANSFORMSTEPSFACTORY_HPP
#define QUICC_TRANSFORM_TRANSFORMSTEPSFACTORY_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Enums/VectorFormulation.hpp"
#include "QuICC/TransformConfigurators/ITransformSteps.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Transform {

   template <typename TSteps> std::shared_ptr<ITransformSteps> tryTransformSteps(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme);

   std::shared_ptr<ITransformSteps> createTransformSteps(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme);

   template <typename TSteps> std::shared_ptr<ITransformSteps> tryTransformSteps(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme)
   {
      if(TSteps::applicable(spScheme))
      {
         return std::make_shared<TSteps>(spScheme);
      } else
      {
         return nullptr;
      }
   }
}
}

#endif // QUICC_TRANSFORM_TRANSFORMSTEPSFACTORY_HPP
