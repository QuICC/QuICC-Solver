/**
 * @file ForwardTransformGrouper.hpp
 * @brief This class defines the forward transform grouping algorithm
 */

#ifndef QUICC_TRANSFORM_FORWARDTRANSFORMGROUPER_HPP
#define QUICC_TRANSFORM_FORWARDTRANSFORMGROUPER_HPP

// Configuration includes
//
#include "QuICC/TypeSelectors/TransformCommSelector.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/PhysicalKernels/IPhysicalKernel.hpp"
#include "QuICC/TransformGroupers/IForwardGrouper.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the forward transform grouping algorithm
    */
   template <typename TConfigurator> class ForwardTransformGrouper: public IForwardGrouper
   {
      public:
         /**
          * @brief Constructor
          */
         ForwardTransformGrouper();

         /**
          * @brief Destructor
          */
         ~ForwardTransformGrouper();

         /**
          * @brief Setup the full forward transform structure for the transform grouping algorithm for the equations
          *
          * @param scalars Vector of scalar fields
          * @param vectors Vector of vector fields
          * @param coord   Transform coord
          */
         virtual void transform(std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalars, std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectors, std::map<std::size_t, Physical::Kernel::SharedIPhysicalKernel>& kernels, TransformCoordinatorType& coord);

         /**
          * @brief Get the number of required buffer packs for the first exchange
          *
          * @param integratorTree Transform integrator tree
          */
         virtual ArrayI packs1D(const std::vector<TransformTree>& integratorTree);

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param integratorTree Transform integrator tree
          */
         virtual ArrayI packs2D(const std::vector<TransformTree>& integratorTree);

      protected:
         /**
          * @brief Setup grouped spectral exchange communication between spectral and first transform
          *
          * @param tree Transform tree describing what fields and what operators to apply
          * @param cood Transform coordinator holding communicators and transforms
          */
         void setupGroupedSpectralCommunication(const TransformTree& tree, TransformCoordinatorType& coord);

         /**
          * @brief Setup grouped first exchange communication between first and second transform
          *
          * @param tree Transform tree describing what fields and what operators to apply
          * @param cood Transform coordinator holding communicators and transforms
          */
         void setupGrouped1DCommunication(TransformCoordinatorType& coord);

         /**
          * @brief Setup grouped second exchange communication between second and third transform
          *
          * @param tree Transform tree describing what fields and what operators to apply
          * @param cood Transform coordinator holding communicators and transforms
          */
         void setupGrouped2DCommunication(TransformCoordinatorType& coord);

         /**
          * @brief Storage for the size of the grouped first exchange communication
          */
         int mGroupedPacks1D;

         /**
          * @brief Storage for the size of the grouped second exchange communication
          */
         int mGroupedPacks2D;

      private:
   };

   template <typename TConfigurator> ForwardTransformGrouper<TConfigurator>::ForwardTransformGrouper()
      : mGroupedPacks1D(-1), mGroupedPacks2D(-1)
   {
      this->split = TConfigurator::SplitLocation;
   }

   template <typename TConfigurator> ForwardTransformGrouper<TConfigurator>::~ForwardTransformGrouper()
   {
   }

   template <typename TConfigurator> inline void ForwardTransformGrouper<TConfigurator>::transform(std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalars, std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectors, std::map<std::size_t, Physical::Kernel::SharedIPhysicalKernel>& kernels, TransformCoordinatorType& coord)
   {
      // Setup the first exchange communication step for scalar fields
      this->setupGrouped2DCommunication(coord);

      //
      // Compute nonlinear interaction
      // ... and first step of forward transform
      //
      for(auto it = coord.forwardTree().cbegin(); it != coord.forwardTree().cend(); ++it)
      {
         if(it->isActive())
         {
            // Get kernel
            auto kernelIt = kernels.find(it->name());

            // Transform scalar equation variable
            if(it->comp<FieldComponents::Physical::Id>() == FieldComponents::Physical::SCALAR)
            {
               auto scalIt = scalars.find(it->name());

               // Compute first step of transform for scalar fields
               TConfigurator::firstStep(*it, scalIt->second, kernelIt->second, coord);

            }
            // Transform vector equation
            else
            {
               auto vectIt = vectors.find(it->name());

               // Compute first step of transform for vector fields
               TConfigurator::firstStep(*it, vectIt->second, kernelIt->second, coord);
            }
         }
      }

      // Initiate the first exchange communication step for scalar fields
      TConfigurator::template initiateCommunication<Dimensions::Transform::TRA3D>(coord);

      // Setup the second exchange communication step for scalar fields
      this->setupGrouped1DCommunication(coord);

      //
      // ... and second step of forward transform
      //
      for(auto it = coord.forwardTree().cbegin(); it != coord.forwardTree().cend(); ++it)
      {
         if(it->isActive())
         {
            // Transform scalar equation variable
            if(it->comp<FieldComponents::Physical::Id>() == FieldComponents::Physical::SCALAR)
            {
               auto scalIt = scalars.find(it->name());

               // Compute second step of transform for scalar fields
               TConfigurator::secondStep(*it, scalIt->second, coord);

            }
            // Transform vector equation
            else
            {
               auto vectIt = vectors.find(it->name());

               // Compute second step of transform for vector fields
               TConfigurator::secondStep(*it, vectIt->second, coord);
            }
         }
      }

      // Initiate the second exchange communication step for vector fields
      TConfigurator::template initiateCommunication<Dimensions::Transform::TRA2D>(coord);

      //
      // ... and last step of forward transform
      //
      for(auto it = coord.forwardTree().cbegin(); it != coord.forwardTree().cend(); ++it)
      {
         if(it->isActive())
         {
            // Transform scalar equation variable
            if(it->comp<FieldComponents::Physical::Id>() == FieldComponents::Physical::SCALAR)
            {
               auto scalIt = scalars.find(it->name());

               // Setup the spectral exchange communication step for scalar fields
               this->setupGroupedSpectralCommunication(*it, coord);

               // Compute last step of transform for scalar fields
               TConfigurator::lastStep(*it, scalIt->second, coord);
               // Initiate the spectral exchange communication step for scalar fields
               TConfigurator::template initiateCommunication<Dimensions::Transform::TRA1D>(coord);

               // Compute spectral step of transform for scalar fields
               TConfigurator::spectralStep(*it, scalIt->second, coord);

            }
            // Transform vector equation
            else
            {
               auto vectIt = vectors.find(it->name());

               // Setup the spectral exchange communication step for vector fields
               this->setupGroupedSpectralCommunication(*it, coord);

               // Compute last step of transform for vector fields
               TConfigurator::lastStep(*it, vectIt->second, coord);
               // Initiate the spectral exchange communication step for vector fields
               TConfigurator::template initiateCommunication<Dimensions::Transform::TRA1D>(coord);

               // Compute spectral step of transform for vector fields
               TConfigurator::spectralStep(*it, vectIt->second, coord);
            }
         }
      }
   }

   template <typename TConfigurator> void ForwardTransformGrouper<TConfigurator>::setupGroupedSpectralCommunication(const TransformTree& tree, TransformCoordinatorType& coord)
   {
      const int packs = this->mNamedPacks1D.at(std::make_pair(tree.name(), tree.comp<FieldComponents::Physical::Id>()));

      TConfigurator::template setupCommunication<Dimensions::Transform::TRA1D>(packs, coord);
   }

   template <typename TConfigurator> void ForwardTransformGrouper<TConfigurator>::setupGrouped1DCommunication(TransformCoordinatorType& coord)
   {
      const int packs = this->mGroupedPacks1D;

      TConfigurator::template setupCommunication<Dimensions::Transform::TRA2D>(packs, coord);
   }

   template <typename TConfigurator> void ForwardTransformGrouper<TConfigurator>::setupGrouped2DCommunication(TransformCoordinatorType& coord)
   {
      const int packs = this->mGroupedPacks2D;

      TConfigurator::template setupCommunication<Dimensions::Transform::TRA3D>(packs, coord);
   }

   template <typename TConfigurator> ArrayI ForwardTransformGrouper<TConfigurator>::packs1D(const std::vector<TransformTree>& integratorTree)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks1D(integratorTree);

      // Store the number of grouped packs
      this->mGroupedPacks1D = packs(0);

      return packs;
   }

   template <typename TConfigurator> ArrayI ForwardTransformGrouper<TConfigurator>::packs2D(const std::vector<TransformTree>& integratorTree)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks2D(integratorTree);

      // Store the number of grouped packs
      this->mGroupedPacks2D = packs(0);

      return packs;
   }

}
}

#endif // QUICC_TRANSFORM_FORWARDTRANSFORMGROUPER_HPP
