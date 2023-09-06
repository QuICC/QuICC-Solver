/**
 * @file BackwardSerialConfigurator.hpp
 * @brief This defines the backward transform serial operations
 */

#ifndef QUICC_TRANSFORM_BACKWARDSERIALCONFIGURATOR_HPP
#define QUICC_TRANSFORM_BACKWARDSERIALCONFIGURATOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/TypeSelectors/TransformCommSelector.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/TransformConfigurators/TransformTree.hpp"
#include "QuICC/TransformConfigurators/BackwardConfigurator.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the backward transform serial operations
    */
   class BackwardSerialConfigurator: public BackwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splitting::Locations::Id  SplitLocation = Splitting::Locations::NONE;

         /**
          * @brief Compute the spectral step in the backward transform
          *
          * @param tree       Transform projector tree
          * @param rVariable  Variable corresponding to the name
          * @param coord      Transform coordinator
          *
          * \tparam TVariable Type of the physical variable
          */
         template <typename TVariable> static void spectralStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord);

         /**
          * @brief Compute the first step in the backward transform
          *
          * @param tree       Transform projector tree
          * @param rVariable  Variable corresponding to the name
          * @param coord      Transform coordinator
          *
          * \tparam TVariable Type of the physical variable
          */
         template <typename TVariable> static void firstStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord);

         /**
          * @brief Compute the second step in the backward transform
          *
          * @param tree       Transform projector tree
          * @param rVariable  Variable corresponding to the name
          * @param coord      Transform coordinator
          *
          * \tparam TVariable Type of the physical variable
          */
         template <typename TVariable> static void secondStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord);

         /**
          * @brief Compute the last step in the backward transform
          *
          * @param tree       Transform projector tree
          * @param rVariable  Variable corresponding to the name
          * @param coord      Transform coordinator
          *
          * \tparam TVariable Type of the physical variable
          */
         template <typename TVariable> static void lastStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord);

         /**
          * @brief Setup exchange communication
          *
          * TId = TRA1D: Communication between Spectral and first transform
          * TId = TRA2D: Communication between first and second transform
          * TId = TRA3D: Communication between second and third transform
          *
          * @param packs Number of components to pack in single communication
          * @param coord Transform coordinator holding communicators and transforms
          *
          * @tparam TId Communication/transpose stage ID
          */
         template <Dimensions::Transform::Id TId> static void setupCommunication(const int packs, TransformCoordinatorType& coord);

         /**
          * @brief Initiate exchange communication
          *
          * TId = TRA1D: Communication between Spectral and first transform
          * TId = TRA2D: Communication between first and second transform
          * TId = TRA3D: Communication between second and third transform
          *
          * @param coord Transform coordinator holding communicators and transforms
          *
          * @tparam TId Communication/transpose stage ID
          */
         template <Dimensions::Transform::Id TId> static void initiateCommunication(TransformCoordinatorType& coord);

      protected:
         /**
          * @brief Empty constructor
          */
         BackwardSerialConfigurator() {};

         /**
          * @brief Empty destructor
          */
         virtual ~BackwardSerialConfigurator() {};

      private:
   };

   template <Dimensions::Transform::Id TId> inline void BackwardSerialConfigurator::setupCommunication(const int, TransformCoordinatorType&)
   {
   }

   template <Dimensions::Transform::Id TId> inline void BackwardSerialConfigurator::initiateCommunication(TransformCoordinatorType&)
   {
   }

   template <typename TVariable> void BackwardSerialConfigurator::spectralStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
   }

   template <typename TVariable> void BackwardSerialConfigurator::firstStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("BwdFirstStep");

      // Iterators for the transforms
      TransformTreeEdge::EdgeType_citerator itSpec;
      TransformTreeEdge::EdgeType_citerator itPhys;

      // Ranges for the vector of edges for the three transforms
      TransformTreeEdge::EdgeType_crange rangeSpec = tree.root().edgeRange();
      TransformTreeEdge::EdgeType_crange rangePhys;

      // Prepare required spectral data
      BackwardConfigurator::prepareSpectral(tree, rVariable, coord);

      if(coord.ss().dimension() == 3)
      {
         // Iterators for the second transforms
         TransformTreeEdge::EdgeType_citerator it2D;

         // Ranges for the vector of edges for the second transforms
         TransformTreeEdge::EdgeType_crange range2D;

         // Loop over first transform
         for(itSpec = rangeSpec.first; itSpec != rangeSpec.second; ++itSpec)
         {
            // Compute first transform
            BackwardConfigurator::project1D(*itSpec, coord);

            range2D = itSpec->edgeRange();
            for(it2D = range2D.first; it2D != range2D.second; ++it2D)
            {
               // Compute second transform
               BackwardConfigurator::project2D(*it2D, coord);

               rangePhys = it2D->edgeRange();
               for(itPhys = rangePhys.first; itPhys != rangePhys.second; ++itPhys)
               {
                  // Prepare physical output data
                  BackwardConfigurator::preparePhysical(tree, *itPhys, rVariable, coord);

                  // Compute third transform
                  BackwardConfigurator::projectND(*itPhys, coord);
               }
            }
         }
      } else if(coord.ss().dimension() == 2)
      {
         // Loop over first transform
         for(itSpec = rangeSpec.first; itSpec != rangeSpec.second; ++itSpec)
         {
            // Compute first transform
            BackwardConfigurator::project1D(*itSpec, coord);

            rangePhys = itSpec->edgeRange();
            for(itPhys = rangePhys.first; itPhys != rangePhys.second; ++itPhys)
            {
               // Prepare physical output data
               BackwardConfigurator::preparePhysical(tree, *itPhys, rVariable, coord);

               // Compute third transform
               BackwardConfigurator::projectND(*itPhys, coord);
            }
         }
      } else if(coord.ss().dimension() == 1)
      {
         throw std::logic_error("1D case is not implemented");
      } else
      {
         throw std::logic_error("Transform with more than 3 dimensions are not implemented");
      }
   }

   template <typename TVariable> void BackwardSerialConfigurator::secondStep(const TransformTree&, TVariable&, TransformCoordinatorType&)
   {
   }

   template <typename TVariable> void BackwardSerialConfigurator::lastStep(const TransformTree&, TVariable&, TransformCoordinatorType&)
   {
   }

}
}

#endif // QUICC_TRANSFORM_BACKWARDSERIALCONFIGURATOR_HPP
