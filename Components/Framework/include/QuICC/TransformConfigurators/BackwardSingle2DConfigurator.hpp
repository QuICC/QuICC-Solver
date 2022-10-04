/**
 * @file BackwardSingle2DConfigurator.hpp
 * @brief This defines the backward transform second exchange single splitting operations
 */

#ifndef QUICC_TRANSFORM_BACKWARDSINGLE2DCONFIGURATOR_HPP
#define QUICC_TRANSFORM_BACKWARDSINGLE2DCONFIGURATOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/TypeSelectors/TransformCommSelector.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/TransformConfigurators/BackwardConfigurator.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the backward transform second exchange single splitting operations
    */
   class BackwardSingle2DConfigurator: public BackwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splitting::Locations::Id  SplitLocation = Splitting::Locations::SECOND;

         /**
          * @brief Compute the first step in the backward transform
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
         BackwardSingle2DConfigurator(){};

         /**
          * @brief Empty destructor
          */
         virtual ~BackwardSingle2DConfigurator(){};

      private:
   };

   template <typename TVariable> void BackwardSingle2DConfigurator::spectralStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("BwdSpectralStep");

      // Prepare required spectral data
      BackwardConfigurator::prepareSpectral(tree, rVariable, coord);
   }

   template <typename TVariable> void BackwardSingle2DConfigurator::firstStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("BwdFirstStep");

      // Iterators for the three transforms
      TransformTreeEdge::EdgeType_citerator it1D;
      TransformTreeEdge::EdgeType_citerator it2D;

      // Ranges for the vector of edges for the three transforms
      TransformTreeEdge::EdgeType_crange range1D = tree.root().edgeRange();
      TransformTreeEdge::EdgeType_crange range2D;

      // Loop over first transform
      for(it1D = range1D.first; it1D != range1D.second; ++it1D)
      {
         // Compute first transform
         BackwardConfigurator::project1D(*it1D, coord);

         range2D = it1D->edgeRange();
         for(it2D = range2D.first; it2D != range2D.second; ++it2D)
         {
            // Compute second transform
            BackwardConfigurator::project2D(*it2D, coord);
         }
      }
   }

   template <typename TVariable> void BackwardSingle2DConfigurator::secondStep(const TransformTree&, TVariable&, TransformCoordinatorType&)
   {
   }

   template <typename TVariable> void BackwardSingle2DConfigurator::lastStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("BwdLastStep");

      // Iterators for the three transforms
      TransformTreeEdge::EdgeType_citerator it1D;
      TransformTreeEdge::EdgeType_citerator it2D;
      TransformTreeEdge::EdgeType_citerator it3D;

      // Ranges for the vector of edges for the three transforms
      TransformTreeEdge::EdgeType_crange range1D = tree.root().edgeRange();
      TransformTreeEdge::EdgeType_crange range2D;
      TransformTreeEdge::EdgeType_crange range3D;

      // Loop over first transform
      for(it1D = range1D.first; it1D != range1D.second; ++it1D)
      {
         range2D = it1D->edgeRange();
         for(it2D = range2D.first; it2D != range2D.second; ++it2D)
         {
            range3D = it2D->edgeRange();
            for(it3D = range3D.first; it3D != range3D.second; ++it3D)
            {
               // Prepare physical output data
               BackwardConfigurator::preparePhysical(tree, *it3D, rVariable, coord);

               // Compute third transform
               BackwardConfigurator::projectND(*it3D, coord);
            }
         }
      }
   }

   template <Dimensions::Transform::Id TId> inline void BackwardSingle2DConfigurator::setupCommunication(const int packs, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("Bwd-setupCommunication");

      if constexpr(TId == Dimensions::Transform::TRA3D)
      {
         coord.communicator().converter<TId>().setupCommunication(packs, TransformDirection::BACKWARD);

         coord.communicator().converter<TId>().prepareBackwardReceive();
      }
   }

   template <Dimensions::Transform::Id TId> inline void BackwardSingle2DConfigurator::initiateCommunication(TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("Bwd-initiateCommunication");

      if constexpr(TId == Dimensions::Transform::TRA3D)
      {
         coord.communicator().converter<TId>().initiateForwardSend();
      }
   }

}
}

#endif // QUICC_TRANSFORM_BACKWARDSINGLE2DCONFIGURATOR_HPP
