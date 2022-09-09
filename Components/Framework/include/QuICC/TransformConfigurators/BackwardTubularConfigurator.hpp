/**
 * @file BackwardTubularConfigurator.hpp
 * @brief This defines the backward transform tubular splitting operations
 */

#ifndef QUICC_TRANSFORM_BACKWARDTUBULARCONFIGURATOR_HPP
#define QUICC_TRANSFORM_BACKWARDTUBULARCONFIGURATOR_HPP

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
#include "QuICC/TransformConfigurators/BackwardConfigurator.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the backward transform tubular splitting operations
    */
   class BackwardTubularConfigurator: public BackwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splitting::Locations::Id  SplitLocation = Splitting::Locations::BOTH;

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
          * @brief Setup first exchange communication
          */
         static void setup1DCommunication(const int packs, TransformCoordinatorType& coord);

         /**
          * @brief Setup second exchange communication
          */
         static void setup2DCommunication(const int packs, TransformCoordinatorType& coord);

         /**
          * @brief Initiate first exchange communication
          */
         static void initiate1DCommunication(TransformCoordinatorType& coord);

         /**
          * @brief Initiate second exchange communication
          */
         static void initiate2DCommunication(TransformCoordinatorType& coord);

      protected:
         /**
          * @brief Empty constructor
          */
         BackwardTubularConfigurator(){};

         /**
          * @brief Empty destructor
          */
         virtual ~BackwardTubularConfigurator(){};

      private:
   };

   template <typename TVariable> void BackwardTubularConfigurator::firstStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("BwdFirstStep");

      // Iterators for the three transforms
      TransformTreeEdge::EdgeType_citerator it1D;

      // Ranges for the vector of edges for the three transforms
      TransformTreeEdge::EdgeType_crange range1D = tree.root().edgeRange();

      // Prepare required spectral data
      BackwardConfigurator::prepareSpectral(tree, rVariable, coord);

      // Loop over first transform
      for(it1D = range1D.first; it1D != range1D.second; ++it1D)
      {
         // Compute first transform
         BackwardConfigurator::project1D(*it1D, coord);
      }
   }

   template <typename TVariable> void BackwardTubularConfigurator::secondStep(const TransformTree& tree, TVariable&, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("BwdSecondStep");

      // Iterators for the three transforms
      TransformTreeEdge::EdgeType_citerator it1D;
      TransformTreeEdge::EdgeType_citerator it2D;

      // Ranges for the vector of edges for the three transforms
      TransformTreeEdge::EdgeType_crange range1D = tree.root().edgeRange();
      TransformTreeEdge::EdgeType_crange range2D;

      // Loop over first transform
      for(it1D = range1D.first; it1D != range1D.second; ++it1D)
      {
         range2D = it1D->edgeRange();
         for(it2D = range2D.first; it2D != range2D.second; ++it2D)
         {
            // Compute second transform
            BackwardConfigurator::project2D(*it2D, coord);
         }
      }
   }

   template <typename TVariable> void BackwardTubularConfigurator::lastStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
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

   inline void BackwardTubularConfigurator::setup1DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("Bwd-setup1DCommunication");

      coord.communicator().converter<Dimensions::Transform::TRA2D>().setupCommunication(packs, TransformDirection::BACKWARD);

      coord.communicator().converter<Dimensions::Transform::TRA2D>().prepareBackwardReceive();
   }

   inline void BackwardTubularConfigurator::setup2DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("Bwd-setup2DCommunication");

      coord.communicator().converter<Dimensions::Transform::TRA3D>().setupCommunication(packs, TransformDirection::BACKWARD);

      coord.communicator().converter<Dimensions::Transform::TRA3D>().prepareBackwardReceive();
   }

   inline void BackwardTubularConfigurator::initiate1DCommunication(TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("Bwd-initiate1DCommunication");

      coord.communicator().converter<Dimensions::Transform::TRA2D>().initiateForwardSend();
   }

   inline void BackwardTubularConfigurator::initiate2DCommunication(TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("Bwd-initiate2DCommunication");

      coord.communicator().converter<Dimensions::Transform::TRA3D>().initiateForwardSend();
   }

}
}

#endif // QUICC_TRANSFORM_BACKWARDTUBULARCONFIGURATOR_HPP
