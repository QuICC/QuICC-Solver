/**
 * @file ForwardTubularConfigurator.hpp
 * @brief This class defines the forward transform tubular splitting operations
 */

#ifndef QUICC_TRANSFORM_FORWARDTUBULARCONFIGURATOR_HPP
#define QUICC_TRANSFORM_FORWARDTUBULARCONFIGURATOR_HPP

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
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the forward transform tubular splitting operations
    */
   class ForwardTubularConfigurator: public ForwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splitting::Locations::Id  SplitLocation = Splitting::Locations::BOTH;

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
         ForwardTubularConfigurator() {};

         /**
          * @brief Empty destructor
          */
         virtual ~ForwardTubularConfigurator() {};

      private:
   };

   inline void ForwardTubularConfigurator::setup1DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      ProfilerMacro_start(Debug::Profiler::FWDTRANSFORM);

      coord.communicator().converter<Dimensions::Transform::TRA2D>().setupCommunication(packs, TransformDirection::FORWARD);

      coord.communicator().converter<Dimensions::Transform::TRA2D>().prepareForwardReceive();

      ProfilerMacro_stop(Debug::Profiler::FWDTRANSFORM);
   }

   inline void ForwardTubularConfigurator::setup2DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      ProfilerMacro_start(Debug::Profiler::FWDTRANSFORM);

      coord.communicator().converter<Dimensions::Transform::TRA3D>().setupCommunication(packs, TransformDirection::FORWARD);

      coord.communicator().converter<Dimensions::Transform::TRA3D>().prepareForwardReceive();

      ProfilerMacro_stop(Debug::Profiler::FWDTRANSFORM);
   }

   inline void ForwardTubularConfigurator::initiate1DCommunication(TransformCoordinatorType& coord)
   {
      ProfilerMacro_start(Debug::Profiler::FWDTRANSFORM);

      coord.communicator().converter<Dimensions::Transform::TRA2D>().initiateBackwardSend();

      ProfilerMacro_stop(Debug::Profiler::FWDTRANSFORM);
   }

   inline void ForwardTubularConfigurator::initiate2DCommunication(TransformCoordinatorType& coord)
   {
      ProfilerMacro_start(Debug::Profiler::FWDTRANSFORM);

      coord.communicator().converter<Dimensions::Transform::TRA3D>().initiateBackwardSend();

      ProfilerMacro_stop(Debug::Profiler::FWDTRANSFORM);
   }

   template <typename TVariable> void ForwardTubularConfigurator::firstStep(const TransformTree& tree, TVariable&, Physical::Kernel::SharedIPhysicalKernel spKernel, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("FwdFirstStep");
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

   template <typename TVariable> void ForwardTubularConfigurator::secondStep(const TransformTree& tree, TVariable&, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("FwdSecondStep");
      ProfilerMacro_start(Debug::Profiler::FWDTRANSFORM);

      // Iterators for the three transforms
      TransformTreeEdge::EdgeType_citerator it2D;
      TransformTreeEdge::EdgeType_citerator it3D;

      // Ranges for the vector of edges for the three transforms
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
         }
      }

      ProfilerMacro_stop(Debug::Profiler::FWDTRANSFORM);
   }

   template <typename TVariable> void ForwardTubularConfigurator::lastStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("FwdLastStep");
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

#endif // QUICC_TRANSFORM_FORWARDTUBULARCONFIGURATOR_HPP
