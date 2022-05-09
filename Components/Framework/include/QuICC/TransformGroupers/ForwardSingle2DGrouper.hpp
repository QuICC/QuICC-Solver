/**
 * @file ForwardSingle2DGrouper.hpp
 * @brief This class defines the forward single grouping exchange grouping algorithm for the second exchange
 */

#ifndef QUICC_TRANSFORM_FORWARDSINGLE2DGROUPER_HPP
#define QUICC_TRANSFORM_FORWARDSINGLE2DGROUPER_HPP

// Configuration includes
//
#include "QuICC/TypeSelectors/TransformCommSelector.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/PhysicalKernels/IPhysicalKernel.hpp"
#include "QuICC/TransformGroupers/IForwardGrouper.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the forward single grouping exchange grouping algorithm for the second exchange
    */
   template <typename TConfigurator> class ForwardSingle2DGrouper : public IForwardGrouper
   {
      public:
         /**
          * @brief Constructor
          */
         ForwardSingle2DGrouper();

         /**
          * @brief Destructor
          */
         ~ForwardSingle2DGrouper();

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
          * @brief Setup grouped first exchange communication
          */
         void setupGrouped1DCommunication(const TransformTree& tree, TransformCoordinatorType& coord);

         /**
          * @brief Setup grouped second exchange communication
          */
         void setupGrouped2DCommunication(TransformCoordinatorType& coord);

         /**
          * @brief Storage for the size of the grouped second exchange communication
          */
         int mGroupedPacks2D;

      private:
   };

   template <typename TConfigurator> ForwardSingle2DGrouper<TConfigurator>::ForwardSingle2DGrouper()
      :mGroupedPacks2D(-1)
   {
      this->split = TConfigurator::SplitLocation;
   }

   template <typename TConfigurator> ForwardSingle2DGrouper<TConfigurator>::~ForwardSingle2DGrouper()
   {
   }

   template <typename TConfigurator> inline void ForwardSingle2DGrouper<TConfigurator>::transform(std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalars, std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectors, std::map<std::size_t, Physical::Kernel::SharedIPhysicalKernel>& kernels, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("transformFwdS2D");

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

               // Transform vector equation
            } else
            {
               auto vectIt = vectors.find(it->name());

               // Compute first step of transform for vector fields
               TConfigurator::firstStep(*it, vectIt->second, kernelIt->second, coord);
            }
         }
      }

      // Initiate the first exchange communication step for vector fields
      TConfigurator::initiate2DCommunication(coord);

      //
      // ... and second and last step of forward transform
      //
      for(auto it = coord.forwardTree().cbegin(); it != coord.forwardTree().cend(); ++it)
      {
         if(it->isActive())
         {
            // Transform scalar equation variable
            if(it->comp<FieldComponents::Physical::Id>() == FieldComponents::Physical::SCALAR)
            {
               auto scalIt = scalars.find(it->name());

               // Setup the second exchange communication step for scalar fields
               this->setupGrouped1DCommunication(*it, coord);

               // Compute second step of transform for scalar fields
               TConfigurator::secondStep(*it, scalIt->second, coord);
               // Initiate the second exchange communication step for scalar fields
               TConfigurator::initiate1DCommunication(coord);

               // Compute last step of transform for scalar fields
               TConfigurator::lastStep(*it, scalIt->second, coord);

               // Transform vector equation
            } else
            {
               auto vectIt = vectors.find(it->name());

               // Setup the second exchange communication step for vector fields
               this->setupGrouped1DCommunication(*it, coord);

               // Compute second step of transform for vector fields
               TConfigurator::secondStep(*it, vectIt->second, coord);
               // Initiate the second exchange communication step for vector fields
               TConfigurator::initiate1DCommunication(coord);

               // Compute last step of transform for vector fields
               TConfigurator::lastStep(*it, vectIt->second, coord);
            }
         }
      }
   }

   template <typename TConfigurator> void ForwardSingle2DGrouper<TConfigurator>::setupGrouped1DCommunication(const TransformTree& tree, TransformCoordinatorType& coord)
   {
      int packs = this->mNamedPacks1D.at(std::make_pair(tree.name(), tree.comp<FieldComponents::Physical::Id>()));

      TConfigurator::setup1DCommunication(packs, coord);
   }

   template <typename TConfigurator> void ForwardSingle2DGrouper<TConfigurator>::setupGrouped2DCommunication(TransformCoordinatorType& coord)
   {
      TConfigurator::setup2DCommunication(this->mGroupedPacks2D, coord);
   }

   template <typename TConfigurator> ArrayI ForwardSingle2DGrouper<TConfigurator>::packs1D(const std::vector<TransformTree>& integratorTree)
   {
      return this->namePacks1D(integratorTree);
   }

   template <typename TConfigurator> ArrayI ForwardSingle2DGrouper<TConfigurator>::packs2D(const std::vector<TransformTree>& integratorTree)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks2D(integratorTree);

      // Store the number of grouped packs
      this->mGroupedPacks2D = packs(0);

      return packs;
   }

}
}

#endif // QUICC_TRANSFORM_FORWARDSINGLE2DGROUPER_HPP
