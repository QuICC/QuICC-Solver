/**
 * @file IBackwardGrouper.hpp
 * @brief This class defines some basic forward transform grouping tools in 2D space
 */

#ifndef QUICC_TRANSFORM_IBACKWARDGROUPER_HPP
#define QUICC_TRANSFORM_IBACKWARDGROUPER_HPP

// Configuration includes
//
#include "QuICC/TypeSelectors/TransformCommSelector.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"

// System includes
//
#include <map>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/TransformConfigurators/TransformTree.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the serial forward transform grouping algorithm in 2D space
    */
   class IBackwardGrouper
   {
      public:
         /// Typedef for field and component ID
         typedef std::pair<std::size_t,FieldComponents::Spectral::Id> FieldIdType;

         /**
          * @brief Setup the full backward transform structure
          *
          * @param scalars Vector of scalar fields
          * @param vectors Vector of vector fields
          * @param coord   Transform coord
          */
         virtual void transform(std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalars, std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectors, TransformCoordinatorType& coord) = 0;

         /**
          * @brief Get the number of required buffer packs for the first exchange
          *
          * @param projectorTree Transform projector tree
          */
         virtual ArrayI packs1D(const std::vector<TransformTree>& projectorTree) = 0;

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param projectorTree Transform projector tree
          */
         virtual ArrayI packs2D(const std::vector<TransformTree>& projectorTree) = 0;

         /**
          * @brief Location of the split in the configurator
          */
         Splitting::Locations::Id split;

      protected:
         /**
          * @brief Get and set the name pack numbers for the first exchange
          *
          * @param projectorTree Transform projector tree
          */
         ArrayI namePacks1D(const std::vector<TransformTree>& projectorTree);

         /**
          * @brief Get the grouped pack number for the first exchange
          *
          * @param projectorTree Transform projector tree
          */
         ArrayI groupPacks1D(const std::vector<TransformTree>& projectorTree);

         /**
          * @brief Get and set the named pack numbers for the second exchange
          *
          * @param projectorTree Transform projector tree
          */
         ArrayI namePacks2D(const std::vector<TransformTree>& projectorTree);

         /**
          * @brief Get the grouped pack number for the second exchange
          *
          * @param projectorTree Transform projector tree
          */
         ArrayI groupPacks2D(const std::vector<TransformTree>& projectorTree);

         /**
          * @brief Storage for named packet sizes for the first exchange
          */
         std::map<FieldIdType, int>  mNamedPacks1D;

         /**
          * @brief Storage for named packet sizes for the second exchange
          */
         std::map<FieldIdType, int>  mNamedPacks2D;

         /**
          * @brief Empty constructor
          */
         IBackwardGrouper();

         /**
          * @brief Empty destructor
          */
         virtual ~IBackwardGrouper();

      private:
   };

   /// Typdef for a smart reference counting pointer to a backward grouper base
   typedef std::shared_ptr<IBackwardGrouper>   SharedIBackwardGrouper;

}
}

#endif // QUICC_TRANSFORM_IBACKWARDGROUPER2D_HPP
