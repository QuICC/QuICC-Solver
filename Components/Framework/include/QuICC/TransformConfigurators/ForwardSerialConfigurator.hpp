/**
 * @file ForwardSerialConfigurator.hpp
 * @brief This class defines the forward transform serial operations
 */

#ifndef QUICC_TRANSFORM_FORWARDSERIALCONFIGURATOR_HPP
#define QUICC_TRANSFORM_FORWARDSERIALCONFIGURATOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/TransformConfigurators/TransformTree.hpp"
#include "QuICC/TransformConfigurators/ForwardConfigurator.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the forward transform serial operations
    */
   class ForwardSerialConfigurator: public ForwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splitting::Locations::Id  SplitLocation = Splitting::Locations::NONE;

         /**
          * @brief First step in transform, including the nonlinear interaction
          *
          * @param spKernel  Shared physical kernel
          * @param rVariable Variable corresponding to the name
          * @param coord     Transform coordinator
          */
         template <typename TVariable> static void firstStep(const TransformTree& tree, TVariable& rVariable, Physical::Kernel::SharedIPhysicalKernel spKernel, TransformCoordinatorType& coord);

         /**
          * @brief Second step in transform
          *
          * @param rVariable Variable corresponding to the name
          * @param coord     Transform coordinator
          */
         template <typename TVariable> static void secondStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord);

         /**
          * @brief Last step in transform
          *
          * @param rVariable Variable corresponding to the name
          * @param coord     Transform coordinator
          */
         template <typename TVariable> static void lastStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord);

         /**
          * @brief Spectral step in transform
          *
          * @param rVariable Variable corresponding to the name
          * @param coord     Transform coordinator
          */
         template <typename TVariable> static void spectralStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord);

         /**
          * @brief exchange communication setup
          *
          * TId = TRA1D: Communication between first and Spectral transform
          * TId = TRA2D: Communication between second and first transform
          * TId = TRA3D: Communication between third and second transform
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
          * TId = TRA1D: Communication between first and Spectral transform
          * TId = TRA2D: Communication between second and first transform
          * TId = TRA3D: Communication between third and second transform
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
         ForwardSerialConfigurator() {};

         /**
          * @brief Empty destructor
          */
         virtual ~ForwardSerialConfigurator() {};

      private:
   };

   template <Dimensions::Transform::Id TId> inline void ForwardSerialConfigurator::setupCommunication(const int, TransformCoordinatorType&)
   {
   }

   template <Dimensions::Transform::Id TId> inline void ForwardSerialConfigurator::initiateCommunication(TransformCoordinatorType&)
   {
   }

   template <typename TVariable> void ForwardSerialConfigurator::firstStep(const TransformTree& tree, TVariable& rVariable, Physical::Kernel::SharedIPhysicalKernel spKernel, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("FwdFirstStep");

      // Iterators for the transforms
      TransformTreeEdge::EdgeType_citerator itSpec;
      TransformTreeEdge::EdgeType_citerator itPhys;

      // Ranges for the vector of edges for the three transforms
      TransformTreeEdge::EdgeType_crange rangeSpec;
      TransformTreeEdge::EdgeType_crange rangePhys = tree.root().edgeRange();

      // Compute the physical space kernel
      ForwardConfigurator::nonlinearTerm(tree, spKernel, coord);

      if(coord.ss().dimension() == 3)
      {
         // Iterators for the second transforms
         TransformTreeEdge::EdgeType_citerator it2D;

         // Ranges for the vector of edges for the second transforms
         TransformTreeEdge::EdgeType_crange range2D;

         // Loop over first transform
         for(itPhys = rangePhys.first; itPhys != rangePhys.second; ++itPhys)
         {
            // Compute third transform
            ForwardConfigurator::integrateND(*itPhys, coord);

            range2D = itPhys->edgeRange();
            for(it2D = range2D.first; it2D != range2D.second; ++it2D)
            {
               // Compute second transform
               ForwardConfigurator::integrate2D(*it2D, coord);

               rangeSpec = it2D->edgeRange();
               for(itSpec = rangeSpec.first; itSpec != rangeSpec.second; ++itSpec)
               {
                  // Compute third transform
                  ForwardConfigurator::integrate1D(*itSpec, coord);

                  // Update equation
                  ForwardConfigurator::updateEquation(*itSpec, rVariable, coord);
               }
            }
         }
      } else if(coord.ss().dimension() == 2)
      {
         // Loop over physical transform
         for(itPhys = rangePhys.first; itPhys != rangePhys.second; ++itPhys)
         {
            // Compute physical transform
            ForwardConfigurator::integrateND(*itPhys, coord);

            rangeSpec = itPhys->edgeRange();
            for(itSpec = rangeSpec.first; itSpec != rangeSpec.second; ++itSpec)
            {
               // Compute third transform
               ForwardConfigurator::integrate1D(*itSpec, coord);

               // Update equation
               ForwardConfigurator::updateEquation(*itSpec, rVariable, coord);
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

   template <typename TVariable> void ForwardSerialConfigurator::secondStep(const TransformTree&, TVariable&, TransformCoordinatorType&)
   {
      // No need for a second step
   }

   template <typename TVariable> void ForwardSerialConfigurator::lastStep(const TransformTree&, TVariable&, TransformCoordinatorType&)
   {
      // No need for a last step
   }

   template <typename TVariable> void ForwardSerialConfigurator::spectralStep(const TransformTree&, TVariable&, TransformCoordinatorType&)
   {
      // No need for a spectral step
   }

}
}

#endif // QUICC_TRANSFORM_FORWARDSERIALCONFIGURATOR_HPP
