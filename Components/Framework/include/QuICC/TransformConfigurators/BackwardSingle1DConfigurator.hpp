/**
 * @file BackwardSingle1DConfigurator.hpp
 * @brief This defines the backward transform first exchange single splitting operations
 */

#ifndef QUICC_TRANSFORM_BACKWARDSINGLE1DCONFIGURATOR_HPP
#define QUICC_TRANSFORM_BACKWARDSINGLE1DCONFIGURATOR_HPP

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
    * @brief This class defines the backward transform first exchange single splitting operations
    */
   class BackwardSingle1DConfigurator: public BackwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splitting::Locations::Id  SplitLocation = Splitting::Locations::FIRST;

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
         BackwardSingle1DConfigurator(){};

         /**
          * @brief Empty destructor
          */
         virtual ~BackwardSingle1DConfigurator(){};

      private:
   };

   template <typename TVariable> void BackwardSingle1DConfigurator::firstStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("BwdFirstStep");

      // Iterators for the three transforms
      TransformTreeEdge::EdgeType_citerator itSpec;

      // Ranges for the vector of edges for the three transforms
      TransformTreeEdge::EdgeType_crange rangeSpec = tree.root().edgeRange();

      // Prepare required spectral data
      BackwardConfigurator::prepareSpectral(tree, rVariable, coord);

      // Loop over first transform
      for(itSpec = rangeSpec.first; itSpec != rangeSpec.second; ++itSpec)
      {
         // Compute first transform
         BackwardConfigurator::project1D(*itSpec, coord);
      }
   }

   template <typename TVariable> void BackwardSingle1DConfigurator::secondStep(const TransformTree&, TVariable&, TransformCoordinatorType&)
   {
   }

   template <typename TVariable> void BackwardSingle1DConfigurator::lastStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("BwdLastStep");

      // Iterators for the transforms
      TransformTreeEdge::EdgeType_citerator itSpec;
      TransformTreeEdge::EdgeType_citerator itPhys;

      // Ranges for the vector of edges for the transforms
      TransformTreeEdge::EdgeType_crange rangeSpec = tree.root().edgeRange();
      TransformTreeEdge::EdgeType_crange rangePhys;

      if(coord.ss().dimension() == 3)
      {
         // Iterators for the second transforms
         TransformTreeEdge::EdgeType_citerator it2D;

         // Ranges for the vector of edges for the second transforms
         TransformTreeEdge::EdgeType_crange range2D;

         // Loop over first transform
         for(itSpec = rangeSpec.first; itSpec != rangeSpec.second; ++itSpec)
         {
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
            rangePhys = itSpec->edgeRange();
            for(itPhys = rangePhys.first; itPhys != rangePhys.second; ++itPhys)
            {
               // Prepare physical output data
               BackwardConfigurator::preparePhysical(tree, *itPhys, rVariable, coord);

               // Compute third transform
               BackwardConfigurator::projectND(*itPhys, coord);
            }
         }
      } else
      {
         throw std::logic_error("Configurator cannot be used with less than 2 dimensions");
      }
   }

   inline void BackwardSingle1DConfigurator::setup1DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("Bwd-setup1DCommunication");

      coord.communicator().converter<Dimensions::Transform::TRA2D>().setupCommunication(packs, TransformDirection::BACKWARD);

      coord.communicator().converter<Dimensions::Transform::TRA2D>().prepareBackwardReceive();
   }

   inline void BackwardSingle1DConfigurator::setup2DCommunication(const int, TransformCoordinatorType&)
   {
   }

   inline void BackwardSingle1DConfigurator::initiate1DCommunication(TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("Bwd-initiate1DCommunication");

      coord.communicator().converter<Dimensions::Transform::TRA2D>().initiateForwardSend();
   }

   inline void BackwardSingle1DConfigurator::initiate2DCommunication(TransformCoordinatorType&)
   {
   }

}
}

#endif // QUICC_TRANSFORM_BACKWARDSINGLE1DCONFIGURATOR_HPP
