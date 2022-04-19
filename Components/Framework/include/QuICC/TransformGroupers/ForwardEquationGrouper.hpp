/**
 * @file ForwardEquationGrouper.hpp
 * @brief This class defines a simple equation wise forward transform grouping algorithm (serial algorithm)
 */
#ifndef QUICC_TRANSFORM_FORWARDEQUATIONGROUPER_HPP
#define QUICC_TRANSFORM_FORWARDEQUATIONGROUPER_HPP

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

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines a simple equation wise forward transform grouping algorithm (serial algorithm)
    */
   template <typename TConfigurator> class ForwardEquationGrouper: public IForwardGrouper
   {
      public:
         /**
          * @brief Constructor
          */
         ForwardEquationGrouper();

         /**
          * @brief Destructor
          */
         ~ForwardEquationGrouper();

         /**
          * @brief Setup the full forward transform structure for the transform grouping algorithm for the equations
          *
          * @param scalars Vector of scalar fields
          * @param vectors Vector of vector fields
          * @param coord   Transform coord
          */
         void transform(std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalars, std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectors, std::map<std::size_t, Physical::Kernel::SharedIPhysicalKernel>& kernels, TransformCoordinatorType& coord) final;

         /**
          * @brief Get the number of required buffer packs for the first exchange
          *
          * @param integratorTree Transform integrator tree
          */
         ArrayI packs1D(const std::vector<TransformTree>& integratorTree) final;

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param integratorTree Transform integrator tree
          */
         ArrayI packs2D(const std::vector<TransformTree>& integratorTree) final;

      protected:
         /**
          * @brief Setup grouped first exchange communication
          */
         void setupGrouped1DCommunication(const TransformTree& tree, TransformCoordinatorType& coord);

         /**
          * @brief Setup grouped second exchange communication
          */
         void setupGrouped2DCommunication(const TransformTree& tree, TransformCoordinatorType& coord);

      private:
   };

   template <typename TConfigurator> ForwardEquationGrouper<TConfigurator>::ForwardEquationGrouper()
   {
      this->split = TConfigurator::SplitLocation;
   }

   template <typename TConfigurator> ForwardEquationGrouper<TConfigurator>::~ForwardEquationGrouper()
   {
   }

   template <typename TConfigurator> inline void ForwardEquationGrouper<TConfigurator>::transform(std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalars, std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vectors, std::map<std::size_t, Physical::Kernel::SharedIPhysicalKernel>& kernels, TransformCoordinatorType& coord)
   {
      //
      // Compute nonlinear interaction
      // ... and forward transform
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
               // Setup the second exchange communication step for scalar fields
               this->setupGrouped1DCommunication(*it, coord);

               // Compute first step of transform for scalar fields
               TConfigurator::firstStep(*it, scalIt->second, kernelIt->second, coord);
               // Initiate the first exchange communication step for scalar fields
               TConfigurator::initiate2DCommunication(coord);

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

               // Setup the first exchange communication step for vector fields
               this->setupGrouped2DCommunication(*it, coord);
               // Setup the second exchange communication step for vector fields
               this->setupGrouped1DCommunication(*it, coord);

               // Compute first step of transform for vector fields
               TConfigurator::firstStep(*it, vectIt->second, kernelIt->second, coord);
               // Initiate the first exchange communication step for vector fields
               TConfigurator::initiate2DCommunication(coord);

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

   template <typename TConfigurator> void ForwardEquationGrouper<TConfigurator>::setupGrouped1DCommunication(const TransformTree& tree, TransformCoordinatorType& coord)
   {
      int packs = this->mNamedPacks1D.at(std::make_pair(tree.name(), tree.comp<FieldComponents::Physical::Id>()));

      TConfigurator::setup1DCommunication(packs, coord);
   }

   template <typename TConfigurator> void ForwardEquationGrouper<TConfigurator>::setupGrouped2DCommunication(const TransformTree& tree, TransformCoordinatorType& coord)
   {
      int packs = this->mNamedPacks2D.at(std::make_pair(tree.name(), tree.comp<FieldComponents::Physical::Id>()));

      TConfigurator::setup2DCommunication(packs, coord);
   }

   template <typename TConfigurator> ArrayI ForwardEquationGrouper<TConfigurator>::packs1D(const std::vector<TransformTree>& integratorTree)
   {
      return this->namePacks1D(integratorTree);
   }

   template <typename TConfigurator> ArrayI ForwardEquationGrouper<TConfigurator>::packs2D(const std::vector<TransformTree>& integratorTree)
   {
      return this->namePacks2D(integratorTree);
   }

}
}

#endif // QUICC_TRANSFORM_FORWARDEQUATIONGROUPER_HPP
