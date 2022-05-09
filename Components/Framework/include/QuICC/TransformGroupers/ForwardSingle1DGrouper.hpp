/**
 * @file ForwardSingle1DGrouper.hpp
 * @brief This class defines the forward single grouping exchange grouping algorithm for the first exchange
 */

#ifndef QUICC_TRANSFORM_FORWARDSINGLE1DGROUPER_HPP
#define QUICC_TRANSFORM_FORWARDSINGLE1DGROUPER_HPP

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
    * @brief This class defines the forward single grouping exchange grouping algorithm for the first exchange
    */
   template <typename TConfigurator> class ForwardSingle1DGrouper: public IForwardGrouper
   {
      public:
         /**
          * @brief Constructor
          */
         ForwardSingle1DGrouper();

         /**
          * @brief Destructor
          */
         ~ForwardSingle1DGrouper();

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
         void setupGrouped1DCommunication(TransformCoordinatorType& coord);

         /**
          * @brief Setup grouped second exchange communication
          */
         void setupGrouped2DCommunication(const TransformTree& tree, TransformCoordinatorType& coord);

         /**
          * @brief Storage for the size of the grouped first echange communication
          */
         int mGroupedPacks1D;

      private:
   };

   template <typename TConfigurator> ForwardSingle1DGrouper<TConfigurator>::ForwardSingle1DGrouper()
      :mGroupedPacks1D(-1)
   {
      this->split = TConfigurator::SplitLocation;
   }

   template <typename TConfigurator> ForwardSingle1DGrouper<TConfigurator>::~ForwardSingle1DGrouper()
   {
   }

   template <typename TConfigurator> inline void ForwardSingle1DGrouper<TConfigurator>::transform(std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalars, std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectors, std::map<std::size_t, Physical::Kernel::SharedIPhysicalKernel>& kernels, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("transformFwdS1D");

      // Setup the second exchange communication step for scalar fields
      this->setupGrouped1DCommunication(coord);

      //
      // Compute nonlinear interaction
      // ... and firts and second forward transform steps
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

               // Setup the first exchange communication step for scalar fields
               this->setupGrouped2DCommunication(*it, coord);

               // Compute first step of transform for scalar fields
               TConfigurator::firstStep(*it, scalIt->second, kernelIt->second, coord);
               // Initiate the first exchange communication step for scalar fields
               TConfigurator::initiate2DCommunication(coord);

               // Compute second step of transform for scalar fields
               TConfigurator::secondStep(*it, scalIt->second, coord);

               // Transform vector equation
            } else
            {
               auto vectIt = vectors.find(it->name());

               // Setup the first exchange communication step for vector fields
               this->setupGrouped2DCommunication(*it, coord);

               // Compute first step of transform for vector fields
               TConfigurator::firstStep(*it, vectIt->second, kernelIt->second, coord);
               // Initiate the first exchange communication step for vector fields
               TConfigurator::initiate2DCommunication(coord);

               // Compute second step of transform for vector fields
               TConfigurator::secondStep(*it, vectIt->second, coord);
            }
         }
      }

      // Initiate the second exchange communication step for scalar fields
      TConfigurator::initiate1DCommunication(coord);

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

               // Compute last step of transform for scalar fields
               TConfigurator::lastStep(*it, scalIt->second, coord);

               // Transform vector equation
            } else
            {
               auto vectIt = vectors.find(it->name());

               // Compute last step of transform for vector fields
               TConfigurator::lastStep(*it, vectIt->second, coord);
            }
         }
      }
   }

   template <typename TConfigurator> void ForwardSingle1DGrouper<TConfigurator>::setupGrouped1DCommunication(TransformCoordinatorType& coord)
   {
      TConfigurator::setup1DCommunication(this->mGroupedPacks1D, coord);
   }

   template <typename TConfigurator> void ForwardSingle1DGrouper<TConfigurator>::setupGrouped2DCommunication(const TransformTree& tree, TransformCoordinatorType& coord)
   {
      int packs = this->mNamedPacks2D.at(std::make_pair(tree.name(), tree.comp<FieldComponents::Physical::Id>()));

      TConfigurator::setup2DCommunication(packs, coord);
   }

   template <typename TConfigurator> ArrayI ForwardSingle1DGrouper<TConfigurator>::packs1D(const std::vector<TransformTree>& integratorTree)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks1D(integratorTree);

      // Store the number of grouped packs
      this->mGroupedPacks1D = packs(0);

      return packs;
   }

   template <typename TConfigurator> ArrayI ForwardSingle1DGrouper<TConfigurator>::packs2D(const std::vector<TransformTree>& integratorTree)
   {
      return this->namePacks2D(integratorTree);
   }

}
}

#endif // QUICC_TRANSFORM_FORWARDSINGLE1DGROUPER_HPP
