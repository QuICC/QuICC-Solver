/**
 * @file ForwardSingle2DConfigurator.hpp
 * @brief This class defines the forward transform single second exchange splitting operations
 */

#ifndef QUICC_TRANSFORM_FORWARDSINGLE2DCONFIGURATOR_HPP
#define QUICC_TRANSFORM_FORWARDSINGLE2DCONFIGURATOR_HPP

// Configuration includes
//
#include "QuICC/Debug/Profiler/ProfilerMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Framework/Selector/ScalarField.hpp"
#include "QuICC/TransformConfigurators/TransformTree.hpp"
#include "QuICC/TransformConfigurators/ForwardConfigurator.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the forward transform single second exchange splitting operations
    */
   class ForwardSingle2DConfigurator: public ForwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splitting::Locations::Id  SplitLocation = Splitting::Locations::SECOND;

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
          * @brief First exchange communication setup
          */
         static void setup1DCommunication(const int packs, TransformCoordinatorType& coord);

         /**
          * @brief Second Exchange communication setup
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
         ForwardSingle2DConfigurator() {};

         /**
          * @brief Empty destructor
          */
         virtual ~ForwardSingle2DConfigurator() {};

      private:
   };

   inline void ForwardSingle2DConfigurator::setup1DCommunication(const int, TransformCoordinatorType&)
   {
   }

   inline void ForwardSingle2DConfigurator::setup2DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      ProfilerMacro_start(Debug::Profiler::FWDTRANSFORM);

      coord.communicator().converter<Dimensions::Transform::TRA3D>().setupCommunication(packs, TransformDirection::FORWARD);

      coord.communicator().converter<Dimensions::Transform::TRA3D>().prepareForwardReceive();

      ProfilerMacro_stop(Debug::Profiler::FWDTRANSFORM);
   }

   inline void ForwardSingle2DConfigurator::initiate1DCommunication(TransformCoordinatorType&)
   {
   }

   inline void ForwardSingle2DConfigurator::initiate2DCommunication(TransformCoordinatorType& coord)
   {
      ProfilerMacro_start(Debug::Profiler::FWDTRANSFORM);

      coord.communicator().converter<Dimensions::Transform::TRA3D>().initiateBackwardSend();

      ProfilerMacro_stop(Debug::Profiler::FWDTRANSFORM);
   }

   template <typename TVariable> void ForwardSingle2DConfigurator::firstStep(const TransformTree& tree, TVariable&, Physical::Kernel::SharedIPhysicalKernel spKernel, TransformCoordinatorType& coord)
   {
      ProfilerMacro_start(Debug::Profiler::FWDTRANSFORM);

      // Iterators for the three transforms
      TransformTreeEdge::EdgeType_citerator it3D;

      // Ranges for the vector of edges for the three transforms
      TransformTreeEdge::EdgeType_crange range3D = tree.root().edgeRange();

      ProfilerMacro_stop(Debug::Profiler::FWDTRANSFORM);

      // Compute the nonlinear interaction
      ForwardConfigurator::nonlinearTerm(tree, spKernel, coord);

      ProfilerMacro_start(Debug::Profiler::FWDTRANSFORM);

      // Loop over first transform
      for(it3D = range3D.first; it3D != range3D.second; ++it3D)
      {
         // Compute third transform
         ForwardConfigurator::integrateND(*it3D, coord);
      }

      ProfilerMacro_stop(Debug::Profiler::FWDTRANSFORM);
   }

   template <typename TVariable> void ForwardSingle2DConfigurator::secondStep(const TransformTree&, TVariable&, TransformCoordinatorType&)
   {
      // No need for a second step
   }

   template <typename TVariable> void ForwardSingle2DConfigurator::lastStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      ProfilerMacro_start(Debug::Profiler::FWDTRANSFORM);

      // Iterators for the three transforms
      TransformTreeEdge::EdgeType_citerator it1D;
      TransformTreeEdge::EdgeType_citerator it2D;
      TransformTreeEdge::EdgeType_citerator it3D;

      // Ranges for the vector of edges for the three transforms
      TransformTreeEdge::EdgeType_crange range1D;
      TransformTreeEdge::EdgeType_crange range2D;
      TransformTreeEdge::EdgeType_crange range3D = tree.root().edgeRange();

      // Loop over first transform
      for(it3D = range3D.first; it3D != range3D.second; ++it3D)
      {
         range2D = it3D->edgeRange();
         for(it2D = range2D.first; it2D != range2D.second; ++it2D)
         {
            // Compute second transform
            ForwardConfigurator::integrate2D(*it2D, coord);

            range1D = it2D->edgeRange();
            for(it1D = range1D.first; it1D != range1D.second; ++it1D)
            {
               // Compute third transform
               ForwardConfigurator::integrate1D(*it1D, coord);

               // Update equation
               ForwardConfigurator::updateEquation(*it1D, rVariable, coord);
            }
         }
      }

      ProfilerMacro_stop(Debug::Profiler::FWDTRANSFORM);
   }

}
}

#endif // QUICC_TRANSFORM_FORWARDSINGLE2DCONFIGURATOR_HPP
