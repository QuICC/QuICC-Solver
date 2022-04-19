/**
 * @file CoreResolution.hpp
 * @brief Definition of a resolution object for a single CPU
 */

#ifndef QUICC_CORERESOLUTION_HPP
#define QUICC_CORERESOLUTION_HPP

// Configuration includes
//

// System includes
//
#include <vector>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Resolutions/TransformResolution.hpp"

namespace QuICC {

   /**
    * @brief Definition of a resolution object for a single CPU
    */
   class CoreResolution
   {
      public:
         /**
          * @brief Constructor
          *
          * @param transformRes Resolution object for the transforms
          */
         explicit CoreResolution(const std::vector<SharedTransformResolution>& transformRes);

         /**
          * @brief Empty Destructor
          */
         ~CoreResolution();

         /**
          * @brief Get transform resolution for corresponding transform
          *
          * @param id ID of the transform
          */
         SharedCTransformResolution dim(const Dimensions::Transform::Id id) const;

         /**
          * @brief Get number of transforms/dimensions
          */
         int nDim() const;

         /**
          * @brief Indexes have been cleared
          */
         bool isCleared() const;

      protected:

      private:
         /**
          * @brief Resolution information for all transforms
          */
         std::vector<SharedTransformResolution>   mTransforms;
   };

   /// Typedef for a shared pointer to a CoreResolution object
   typedef std::shared_ptr<CoreResolution>   SharedCoreResolution;

   /// Typedef for a shared pointer to a const CoreResolution object
   typedef std::shared_ptr<const CoreResolution>   SharedCCoreResolution;

}

#endif // QUICC_CORERESOLUTION_HPP
