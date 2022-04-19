/**
 * @file TransformStepsFactory.cpp
 * @brief Source of the factory to create transform steps
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/TransformConfigurators/ITransformSteps.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

   ITransformSteps::ITransformSteps(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme)
      : mspScheme(spScheme)
   {
   }

   ITransformSteps::~ITransformSteps()
   {
   }

   const SpatialScheme::ISpatialScheme& ITransformSteps::ss() const
   {
      return *this->mspScheme;
   }

}
}
