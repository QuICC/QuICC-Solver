/**
 * @file BackwardSingle1DGrouper.hpp
 * @brief This class defines the backward single grouping exchange grouping algorithm for the first exchange
 */

#ifndef QUICC_TRANSFORM_BACKWARDSINGLE1DGROUPER_HPP
#define QUICC_TRANSFORM_BACKWARDSINGLE1DGROUPER_HPP

// Configuration includes
//
#include "QuICC/Enums/Dimensions.hpp"
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
    * @brief This class defines the backward single grouping exchange grouping algorithm for the first exchange
    */
   template <typename TConfigurator> class BackwardSingle1DGrouper: public IBackwardGrouper
   {
      public:
         /**
          * @brief Constructor
          */
         BackwardSingle1DGrouper();

         /**
          * @brief Destructor
          */
         ~BackwardSingle1DGrouper();

         /**
          * @brief Setup the full backward transform structure for the first exchange grouping algorithm
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
          * @brief Get the number of required buffer packs for the second exchange
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
         void setupGrouped2DCommunication(const TransformTree& tree, TransformCoordinatorType& coord);

         /**
          * @brief Storage for the size of the grouped first exchange communication
          */
         int mGroupedPacks1D;

      private:
   };

   template <typename TConfigurator> BackwardSingle1DGrouper<TConfigurator>::BackwardSingle1DGrouper()
      : mGroupedPacks1D(-1)
   {
      this->split = TConfigurator::SplitLocation;
   }

   template <typename TConfigurator> BackwardSingle1DGrouper<TConfigurator>::~BackwardSingle1DGrouper()
   {
   }

   template <typename TConfigurator> inline void BackwardSingle1DGrouper<TConfigurator>::transform(std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalars, std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectors, TransformCoordinatorType& coord)
   {
      //
      // Compute backward transform
      //

      // Setup the grouped first exchange communication
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

      // Initiate the first exchange communication step for vector fields
      TConfigurator::template initiateCommunication<Dimensions::Transform::TRA2D>(coord);

      //
      // Compute intermediate and last steps
      //
      for(auto it = coord.backwardTree().cbegin(); it != coord.backwardTree().cend(); ++it)
      {
         if(it->isActive())
         {
            // Transform scalar variable
            if(it->comp<FieldComponents::Spectral::Id>() == FieldComponents::Spectral::SCALAR)
            {
               auto scalIt = scalars.find(it->name());

               // Setup the second exchange communication step for scalar fields
               this->setupGrouped2DCommunication(*it, coord);

               // Compute second step of transform for scalar fields
               TConfigurator::secondStep(*it, scalIt->second, coord);
               // Initiate the second exchange communication step for scalar fields
               TConfigurator::template initiateCommunication<Dimensions::Transform::TRA3D>(coord);

               // Compute last step of transform for scalar fields
               TConfigurator::lastStep(*it, scalIt->second, coord);
            }
            // Transform vector variable
            else
            {
               auto vectIt = vectors.find(it->name());

               // Setup the second exchange communication step for vector fields
               this->setupGrouped2DCommunication(*it, coord);

               // Compute second step of transform for vector fields
               TConfigurator::secondStep(*it, vectIt->second, coord);
               // Initiate the second exchange communication step for vector fields
               TConfigurator::template initiateCommunication<Dimensions::Transform::TRA3D>(coord);

               // Compute last step of transform for vector fields
               TConfigurator::lastStep(*it, vectIt->second, coord);
            }
         }
      }
   }

   template <typename TConfigurator> void BackwardSingle1DGrouper<TConfigurator>::setupGroupedSpectralCommunication(const TransformTree& tree, TransformCoordinatorType& coord)
   {
      const int packs = 1;

      TConfigurator::template setupCommunication<Dimensions::Transform::TRA1D>(packs, coord);
   }

   template <typename TConfigurator> void BackwardSingle1DGrouper<TConfigurator>::setupGrouped1DCommunication(TransformCoordinatorType& coord)
   {
      if(this->mGroupedPacks1D > 0)
      {
         TConfigurator::template setupCommunication<Dimensions::Transform::TRA2D>(this->mGroupedPacks1D, coord);
      }
   }

   template <typename TConfigurator> void BackwardSingle1DGrouper<TConfigurator>::setupGrouped2DCommunication(const TransformTree& tree, TransformCoordinatorType& coord)
   {
      std::pair<std::size_t,FieldComponents::Spectral::Id> id = std::make_pair(tree.name(), tree.comp<FieldComponents::Spectral::Id>());
      if(this->mNamedPacks2D.count(id) == 1)
      {
         TConfigurator::template setupCommunication<Dimensions::Transform::TRA3D>(this->mNamedPacks2D.at(id), coord);
      }
   }

   template <typename TConfigurator> ArrayI BackwardSingle1DGrouper<TConfigurator>::packs1D(const std::vector<TransformTree>& projectorTree)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks1D(projectorTree);

      // Store the number of grouped packs
      this->mGroupedPacks1D = packs(0);

      return packs;
   }

   template <typename TConfigurator> ArrayI BackwardSingle1DGrouper<TConfigurator>::packs2D(const std::vector<TransformTree>& projectorTree)
   {
      return this->namePacks2D(projectorTree);
   }

}
}

#endif // QUICC_TRANSFORM_BACKWARDSINGLE1DGROUPER_HPP
