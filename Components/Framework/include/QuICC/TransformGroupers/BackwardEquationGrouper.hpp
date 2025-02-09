/**
 * @file BackwardEquationGrouper.hpp
 * @brief This class defines a simple equation wise backward transform grouping algorithm (serial algorithm)
 */
#ifndef QUICC_TRANSFORM_BACKWARDDEQUATIONGROUPER_HPP
#define QUICC_TRANSFORM_BACKWARDDEQUATIONGROUPER_HPP

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
#include "QuICC/Debug/DebuggerMacro.h"
#ifdef QUICC_DEBUG
   #include "QuICC/PhysicalNames/Coordinator.hpp"
   #include "QuICC/Tools/IdToHuman.hpp"
#endif // QUICC_DEBUG
#include "QuICC/TransformGroupers/IBackwardGrouper.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines a simple equation wise backward transform grouping algorithm (serial algorithm)
    */
   template <typename TConfigurator> class BackwardEquationGrouper: public IBackwardGrouper
   {
      public:
         /**
          * @brief Constructor
          */
         BackwardEquationGrouper();

         /**
          * @brief Destructor
          */
         ~BackwardEquationGrouper();

         /**
          * @brief Setup the full backward transform structure
          *
          * @param scalars Vector of scalar fields
          * @param vectors Vector of vector fields
          * @param coord   Transform coord
          */
         void transform(std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalars, std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectors, TransformCoordinatorType& coord) final;

         /**
          * @brief Get the number of required buffer packs for the first exchange
          *
          * @param projectorTree Transform projector tree
          */
         ArrayI packs1D(const std::vector<TransformTree>& projectorTree) final;

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param projectorTree Transform projector tree
          */
         ArrayI packs2D(const std::vector<TransformTree>& projectorTree) final;

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
         void setupGrouped1DCommunication(const TransformTree& tree, TransformCoordinatorType& coord);

         /**
          * @brief Setup grouped second exchange communication between third and second transform
          *
          * @param tree Transform tree describing what fields and what operators to apply
          * @param cood Transform coordinator holding communicators and transforms
          */
         void setupGrouped2DCommunication(const TransformTree& tree, TransformCoordinatorType& coord);

      private:
   };

   template <typename TConfigurator> BackwardEquationGrouper<TConfigurator>::BackwardEquationGrouper()
   {
      this->split = TConfigurator::SplitLocation;
   }

   template <typename TConfigurator> BackwardEquationGrouper<TConfigurator>::~BackwardEquationGrouper()
   {
   }

   template <typename TConfigurator> inline void BackwardEquationGrouper<TConfigurator>::transform(std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalars, std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectors, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("transformBwd");
      //
      // Compute backward transform
      //
      for(auto it = coord.backwardTree().cbegin(); it != coord.backwardTree().cend(); ++it)
      {
         if(it->isActive())
         {
            // Transform scalar variable
            if(it->comp<FieldComponents::Spectral::Id>() == FieldComponents::Spectral::SCALAR)
            {
               DebuggerMacro_msg("Backward transform for " + PhysicalNames::Coordinator::tag(it->name()), 5);
               Profiler::RegionFixture<1> fix("transformBwdScalar");

               auto scalIt = scalars.find(it->name());

               // Setup the first exchange communication step for scalar fields
               this->setupGroupedSpectralCommunication(*it, coord);
               // Setup the first exchange communication step for scalar fields
               this->setupGrouped1DCommunication(*it, coord);
               // Setup the second exchange communication step for scalar fields
               this->setupGrouped2DCommunication(*it, coord);

               // Compute spectral step of transform for scalar fields
               TConfigurator::spectralStep(*it, scalIt->second, coord);
               // Initiate the first exchange communication step for scalar fields
               TConfigurator::template initiateCommunication<Dimensions::Transform::TRA1D>(coord);

               // Compute first step of transform for scalar fields
               TConfigurator::firstStep(*it, scalIt->second, coord);
               // Initiate the first exchange communication step for scalar fields
               TConfigurator::template initiateCommunication<Dimensions::Transform::TRA2D>(coord);

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
               DebuggerMacro_msg("Backward transform for " + PhysicalNames::Coordinator::tag(it->name()) + "(" + Tools::IdToHuman::toString(it->comp<FieldComponents::Spectral::Id>()) + ")", 5);
               Profiler::RegionFixture<1> fix("transformBwdVector");

               auto vectIt = vectors.find(it->name());

               // Setup the first exchange communication step for scalar fields
               this->setupGroupedSpectralCommunication(*it, coord);
               // Setup the first exchange communication step for vector fields
               this->setupGrouped1DCommunication(*it, coord);
               // Setup the second exchange communication step for vector fields
               this->setupGrouped2DCommunication(*it, coord);

               // Compute spectral step of transform for vector fields
               TConfigurator::spectralStep(*it, vectIt->second, coord);
               // Initiate the spectral exchange communication step for vector fields
               TConfigurator::template initiateCommunication<Dimensions::Transform::TRA1D>(coord);

               // Compute first step of transform for vector fields
               TConfigurator::firstStep(*it, vectIt->second, coord);
               // Initiate the first exchange communication step for vector fields
               TConfigurator::template initiateCommunication<Dimensions::Transform::TRA2D>(coord);

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

   template <typename TConfigurator> void BackwardEquationGrouper<TConfigurator>::setupGroupedSpectralCommunication(const TransformTree& tree, TransformCoordinatorType& coord)
   {
      const int packs = 1;

      TConfigurator::template setupCommunication<Dimensions::Transform::TRA1D>(packs, coord);
   }

   template <typename TConfigurator> void BackwardEquationGrouper<TConfigurator>::setupGrouped1DCommunication(const TransformTree& tree, TransformCoordinatorType& coord)
   {
      const int packs = this->mNamedPacks1D.at(std::make_pair(tree.name(), tree.comp<FieldComponents::Spectral::Id>()));

      TConfigurator::template setupCommunication<Dimensions::Transform::TRA2D>(packs, coord);
   }

   template <typename TConfigurator> void BackwardEquationGrouper<TConfigurator>::setupGrouped2DCommunication(const TransformTree& tree, TransformCoordinatorType& coord)
   {
      const int packs = this->mNamedPacks2D.at(std::make_pair(tree.name(), tree.comp<FieldComponents::Spectral::Id>()));

      TConfigurator::template setupCommunication<Dimensions::Transform::TRA3D>(packs, coord);
   }

   template <typename TConfigurator> ArrayI BackwardEquationGrouper<TConfigurator>::packs1D(const std::vector<TransformTree>& projectorTree)
   {
      return this->namePacks1D(projectorTree);
   }

   template <typename TConfigurator> ArrayI BackwardEquationGrouper<TConfigurator>::packs2D(const std::vector<TransformTree>& projectorTree)
   {
      return this->namePacks2D(projectorTree);
   }

}
}

#endif // QUICC_TRANSFORM_BACKWARDDEQUATIONGROUPER_HPP
