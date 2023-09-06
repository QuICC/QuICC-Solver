/**
 * @file BackwardTransformGrouper.hpp
 * @brief This class defines the backward transform grouping algorithm
 */

#ifndef QUICC_TRANSFORM_BACKWARDTRANSFORMGROUPER_HPP
#define QUICC_TRANSFORM_BACKWARDTRANSFORMGROUPER_HPP

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
#include "QuICC/TransformGroupers/IBackwardGrouper.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the backward transform grouping algorithm
    */
   template <typename TConfigurator> class BackwardTransformGrouper: public IBackwardGrouper
   {
      public:
         /**
          * @brief Constructor
          */
         BackwardTransformGrouper();

         /**
          * @brief Destructor
          */
         ~BackwardTransformGrouper();

         /**
          * @brief Setup the full backward transform structure for the transform grouping algorithm
          *
          * @param scalars Vector of scalar fields
          * @param vectors Vector of vector fields
          * @param coord   Transform coord
          */
         virtual void transform(std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalars, std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectors, TransformCoordinatorType& coord);

         /**
          * @brief Get the number of required buffer packs for the first exchange
          *
          * @param projectorTree Transform projector tree
          */
         virtual ArrayI packs1D(const std::vector<TransformTree>& projectorTree);

         /**
          * @brief Get the number of required packs for the second exchange
          *
          * @param projectorTree Transform projector tree
          */
         virtual ArrayI packs2D(const std::vector<TransformTree>& projectorTree);

      protected:
         /**
          * @brief Setup grouped spectral exchange communication between first transform and spectra
          *
          * @param tree Transform tree describing what fields and what operators to apply
          * @param cood Transform coordinator holding communicators and transforms
          */
         void setupGroupedSpectralCommunication(const TransformTree& tree, TransformCoordinatorType& coord);

         /**
          * @brief Setup grouped first exchange communication between second and first transform
          *
          * @param tree Transform tree describing what fields and what operators to apply
          * @param cood Transform coordinator holding communicators and transforms
          */
         void setupGrouped1DCommunication(TransformCoordinatorType& coord);

         /**
          * @brief Setup grouped second exchange communication between third and second transform
          *
          * @param tree Transform tree describing what fields and what operators to apply
          * @param cood Transform coordinator holding communicators and transforms
          */
         void setupGrouped2DCommunication(TransformCoordinatorType& coord);

         /**
          * @brief Storage for the size of the grouped communication for the first exchange
          */
         int mGroupedPacks1D;

         /**
          * @brief Storage for the size of the grouped communication for the second exchange
          */
         int mGroupedPacks2D;

      private:
   };

   template <typename TConfigurator> BackwardTransformGrouper<TConfigurator>::BackwardTransformGrouper()
      : mGroupedPacks1D(-1), mGroupedPacks2D(-1)
   {
      this->split = TConfigurator::SplitLocation;
   }

   template <typename TConfigurator> BackwardTransformGrouper<TConfigurator>::~BackwardTransformGrouper()
   {
   }

   template <typename TConfigurator> inline void BackwardTransformGrouper<TConfigurator>::transform(std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalars, std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectors, TransformCoordinatorType& coord)
   {
      //
      // Compute backward transform
      //

      // Setup the 1D grouped first exchange communication
      this->setupGrouped1DCommunication(coord);

      //
      // Compute first step of backward transform
      //
      for(auto it = coord.backwardTree().cbegin(); it != coord.backwardTree().cend(); ++it)
      {
         if(it->isActive())
         {
            // Transform scalar variable
            if(it->comp<FieldComponents::Spectral::Id>() == FieldComponents::Spectral::SCALAR)
            {
               auto scalIt = scalars.find(it->name());

               // Setup the first exchange communication step for scalar fields
               this->setupGroupedSpectralCommunication(*it, coord);

               // Compute spectral step of transform for scalar fields
               TConfigurator::spectralStep(*it, scalIt->second, coord);
               // Initiate the first exchange communication step for scalar fields
               TConfigurator::template initiateCommunication<Dimensions::Transform::TRA1D>(coord);

               // Compute first step of transform for scalar fields
               TConfigurator::firstStep(*it, scalIt->second, coord);
            }
            // Transform vector variable
            else
            {
               auto vectIt = vectors.find(it->name());

               // Setup the first exchange communication step for scalar fields
               this->setupGroupedSpectralCommunication(*it, coord);

               // Compute spectral step of transform for vector fields
               TConfigurator::spectralStep(*it, vectIt->second, coord);
               // Initiate the spectral exchange communication step for vector fields
               TConfigurator::template initiateCommunication<Dimensions::Transform::TRA1D>(coord);

               // Compute first step of transform for vector fields
               TConfigurator::firstStep(*it, vectIt->second, coord);
            }
         }
      }

      // Initiate the grouped first exchange communication
      TConfigurator::template initiateCommunication<Dimensions::Transform::TRA2D>(coord);

      // Setup the grouped second exchange communication
      this->setupGrouped2DCommunication(coord);

      //
      // Compute intermediate step
      //
      for(auto it = coord.backwardTree().cbegin(); it != coord.backwardTree().cend(); ++it)
      {
         if(it->isActive())
         {
            // Transform scalar variable
            if(it->comp<FieldComponents::Spectral::Id>() == FieldComponents::Spectral::SCALAR)
            {
               auto scalIt = scalars.find(it->name());

               // Compute second step of transform for scalar fields
               TConfigurator::secondStep(*it, scalIt->second, coord);
            }
            // Transform vector variable
            else
            {
               auto vectIt = vectors.find(it->name());

               // Compute second step of transform for vector fields
               TConfigurator::secondStep(*it, vectIt->second, coord);
            }
         }
      }

      // Initiate the grouped second exchange communication
      TConfigurator::template initiateCommunication<Dimensions::Transform::TRA3D>(coord);

      //
      // Compute last step
      //
      for(auto it = coord.backwardTree().cbegin(); it != coord.backwardTree().cend(); ++it)
      {
         if(it->isActive())
         {
            // Transform scalar variable
            if(it->comp<FieldComponents::Spectral::Id>() == FieldComponents::Spectral::SCALAR)
            {
               auto scalIt = scalars.find(it->name());

               // Compute last step of transform for scalar fields
               TConfigurator::lastStep(*it, scalIt->second, coord);
            }
            // Transform vector variable
            else
            {
               auto vectIt = vectors.find(it->name());

               // Compute last step of transform for vector fields
               TConfigurator::lastStep(*it, vectIt->second, coord);
            }
         }
      }
   }

   template <typename TConfigurator> void BackwardTransformGrouper<TConfigurator>::setupGroupedSpectralCommunication(const TransformTree& tree, TransformCoordinatorType& coord)
   {
      const int packs = 1;

      TConfigurator::template setupCommunication<Dimensions::Transform::TRA1D>(packs, coord);
   }

   template <typename TConfigurator> void BackwardTransformGrouper<TConfigurator>::setupGrouped1DCommunication(TransformCoordinatorType& coord)
   {
      if(this->mGroupedPacks1D > 0)
      {
         TConfigurator::template setupCommunication<Dimensions::Transform::TRA2D>(this->mGroupedPacks1D, coord);
      }
   }

   template <typename TConfigurator> void BackwardTransformGrouper<TConfigurator>::setupGrouped2DCommunication(TransformCoordinatorType& coord)
   {
      if(this->mGroupedPacks2D > 0)
      {
         TConfigurator::template setupCommunication<Dimensions::Transform::TRA3D>(this->mGroupedPacks2D, coord);
      }
   }

   template <typename TConfigurator> ArrayI BackwardTransformGrouper<TConfigurator>::packs1D(const std::vector<TransformTree>& projectorTree)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks1D(projectorTree);

      // Store the number of grouped packs
      this->mGroupedPacks1D = packs(0);

      return packs;
   }

   template <typename TConfigurator> ArrayI BackwardTransformGrouper<TConfigurator>::packs2D(const std::vector<TransformTree>& projectorTree)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks2D(projectorTree);

      // Store the number of grouped packs
      this->mGroupedPacks2D = packs(0);

      return packs;
   }

}
}

#endif // QUICC_TRANSFORM_BACKWARDTRANSFORMGROUPER_HPP
