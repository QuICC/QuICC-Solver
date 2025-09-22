/**
 * @file TransformStepsFactory.cpp
 * @brief Source of the factory to create transform steps
 */

// System includes
//

// Project includes
//
#include "QuICC/TransformConfigurators/ITransformSteps.hpp"

namespace QuICC {

namespace Transform {

   ITransformSteps::ITransformSteps(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme)
      : mspScheme(spScheme)
   {
   }

   const SpatialScheme::ISpatialScheme& ITransformSteps::ss() const
   {
      return *this->mspScheme;
   }

}
}
