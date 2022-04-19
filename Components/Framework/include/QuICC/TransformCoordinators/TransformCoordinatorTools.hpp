/**
 * @file TransformCoordinatorTools.hpp
 * @brief Implementation of transform coordinator tools
 */

#ifndef QUICC_TRANSFORM_TRANSFORMCOORDINATORTOOLS_HPP
#define QUICC_TRANSFORM_TRANSFORMCOORDINATORTOOLS_HPP

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/TypeSelectors/TransformCommSelector.hpp"
#include "QuICC/TransformConfigurators/TransformTree.hpp"
#include "QuICC/TransformGroupers/IForwardGrouper.hpp"
#include "QuICC/TransformGroupers/IBackwardGrouper.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Timers/StageTimer.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Implementation of transform coordinator tools
    */
   class TransformCoordinatorTools
   {
      public:
         /**
          * @brief Initialise the transform coordinator
          *
          * @param rCoord           Transform coordinator
          * @param spFwdGrouper     Forward transform communication grouper
          * @param spBwdGrouper     Backward transform communication grouper
          * @param forwardTree   Transform integrator tree
          * @param backwardTree    Transform projector tree
          * @param spRes            Shared resolution
          * @param runOptions       Available run options map
          */
         static void init(TransformCoordinatorType& rCoord, SharedIForwardGrouper spFwdGrouper, SharedIBackwardGrouper spBwdGrouper, const std::vector<Transform::TransformTree>& forwardTree, const std::vector<Transform::TransformTree>& backwardTree, SharedResolution spRes, const std::map<std::size_t,NonDimensional::SharedINumber>& runOptions);

      protected:

      private:
         /**
          * @brief Constructor
          */
         TransformCoordinatorTools();

         /**
          * @brief Destructor
          */
         ~TransformCoordinatorTools();
   };

}
}

#endif // QUICC_TRANSFORM_TRANSFORMCOORDINATORTOOLS_HPP
