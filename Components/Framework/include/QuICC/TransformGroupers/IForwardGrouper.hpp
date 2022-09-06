/**
 * @file IForwardGrouper.hpp
 * @brief This class defines some basic forward transform grouping tools in xD Space
 */

#ifndef QUICC_TRANSFORM_IFORWARDGROUPER_HPP
#define QUICC_TRANSFORM_IFORWARDGROUPER_HPP

// Configuration includes
//
#include "QuICC/TypeSelectors/TransformCommSelector.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/PhysicalKernels/IPhysicalKernel.hpp"
#include "QuICC/TransformConfigurators/TransformTree.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the serial forward transform grouping algorithm in xD space
    */
   class IForwardGrouper
   {
      public:
         /// Typedef for field and component ID
         typedef std::pair<std::size_t,FieldComponents::Physical::Id> FieldIdType;

         /**
          * @brief Setup the full forward transform structure for the transform grouping algorithm for the equations
          *
          * @param scalars Vector of scalar fields
          * @param vectors Vector of vector fields
          * @param coord   Transform coord
          */
         virtual void transform(std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalars, std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectors, std::map<std::size_t, Physical::Kernel::SharedIPhysicalKernel>& kernels, TransformCoordinatorType& coord) = 0;

         /**
          * @brief Get the number of required buffer packs for the first exchange
          *
          * @param integratorTree Transform integrator tree
          */
         virtual ArrayI packs1D(const std::vector<TransformTree>& integratorTree) = 0;

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param integratorTree Transform integrator tree
          */
         virtual ArrayI packs2D(const std::vector<TransformTree>& integratorTree) = 0;

         /**
          * @brief Location of the split in the configurator
          */
         Splitting::Locations::Id split;

      protected:
         /**
          * @brief Get and set the name pack numbers for the first exchange
          *
          * @param integratorTree Transform integrator tree
          */
         ArrayI namePacks1D(const std::vector<TransformTree>& integratorTree);

         /**
          * @brief Get the grouped pack number for the first exchange
          *
          * @param integratorTree Transform integrator tree
          */
         ArrayI groupPacks1D(const std::vector<TransformTree>& integratorTree);

         /**
          * @brief Get and set the named pack numbers for the second exchange
          *
          * @param integratorTree Transform integrator tree
          */
         ArrayI namePacks2D(const std::vector<TransformTree>& integratorTree);

         /**
          * @brief Get the grouped pack number for the second exchange
          *
          * @param integratorTree Transform integrator tree
          */
         ArrayI groupPacks2D(const std::vector<TransformTree>& integratorTree);

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
         IForwardGrouper();

         /**
          * @brief Empty destructor
          */
         virtual ~IForwardGrouper();

      private:
   };

   /// Typdef for a smart reference counting pointer to a backward grouper base
   typedef std::shared_ptr<IForwardGrouper>   SharedIForwardGrouper;

}
}

#endif // QUICC_TRANSFORM_IFORWARDGROUPER_HPP
