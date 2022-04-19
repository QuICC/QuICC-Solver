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
#include "QuICC/TransformConfigurators/TransformStepsFactory.hpp"

// Project includes
//
#include "QuICC/TransformConfigurators/ITransformSteps.hpp"
#include "QuICC/TransformConfigurators/CartesianTransformSteps.hpp"
#include "QuICC/TransformConfigurators/AnnulusTransformSteps.hpp"
#include "QuICC/TransformConfigurators/CylinderTransformSteps.hpp"
#include "QuICC/TransformConfigurators/SphereTransformSteps.hpp"
#include "QuICC/TransformConfigurators/ShellTransformSteps.hpp"

namespace QuICC {

namespace Transform {

   std::shared_ptr<ITransformSteps> createTransformSteps(std::shared_ptr<const SpatialScheme::ISpatialScheme> spScheme)
   {
      std::shared_ptr<ITransformSteps> spSteps;

      if(!spSteps) spSteps = tryTransformSteps<SphereTransformSteps>(spScheme);
      if(!spSteps) spSteps = tryTransformSteps<ShellTransformSteps>(spScheme);
      if(!spSteps) spSteps = tryTransformSteps<CartesianTransformSteps>(spScheme);
      if(!spSteps) spSteps = tryTransformSteps<CylinderTransformSteps>(spScheme);
      if(!spSteps) spSteps = tryTransformSteps<AnnulusTransformSteps>(spScheme);

      return spSteps;
   }

}
}
