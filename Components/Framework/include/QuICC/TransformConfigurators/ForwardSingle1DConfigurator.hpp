/**
 * @file ForwardSingle1DConfigurator.hpp
 * @brief This class defines the forward transform single first exchange splitting operations
 */

#ifndef QUICC_TRANSFORM_FORWARDSINGLE1DCONFIGURATOR_HPP
#define QUICC_TRANSFORM_FORWARDSINGLE1DCONFIGURATOR_HPP

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
    * @brief This class defines the forward transform single first exchange splitting operations
    */
   class ForwardSingle1DConfigurator: public ForwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splitting::Locations::Id  SplitLocation = Splitting::Locations::FIRST;

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
         ForwardSingle1DConfigurator() {};

         /**
          * @brief Empty destructor
          */
         virtual ~ForwardSingle1DConfigurator() {};

      private:
   };

   inline void ForwardSingle1DConfigurator::setup1DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("Fwd-setup1DCommunication");

      coord.communicator().converter<Dimensions::Transform::TRA2D>().setupCommunication(packs, TransformDirection::FORWARD);

      coord.communicator().converter<Dimensions::Transform::TRA2D>().prepareForwardReceive();
   }

   inline void ForwardSingle1DConfigurator::setup2DCommunication(const int, TransformCoordinatorType&)
   {
   }

   inline void ForwardSingle1DConfigurator::initiate1DCommunication(TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<2> fix("Fwd-initiate1DCommunication");

      coord.communicator().converter<Dimensions::Transform::TRA2D>().initiateBackwardSend();
   }

   inline void ForwardSingle1DConfigurator::initiate2DCommunication(TransformCoordinatorType&)
   {
   }

   template <typename TVariable> void ForwardSingle1DConfigurator::firstStep(const TransformTree& tree, TVariable&, Physical::Kernel::SharedIPhysicalKernel spKernel, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("FwdFirstStep");

      // Iterators for the transforms
      TransformTreeEdge::EdgeType_citerator it3D;

      // Ranges for the vector of edges for the transforms
      TransformTreeEdge::EdgeType_crange range3D = tree.root().edgeRange();

      // Compute the nonlinear interaction
      ForwardConfigurator::nonlinearTerm(tree, spKernel, coord);

      if(coord.ss().dimension() == 3)
      {
         // Iterators for the second transforms
         TransformTreeEdge::EdgeType_citerator it2D;

         // Ranges for the vector of edges for the second transforms
         TransformTreeEdge::EdgeType_crange range2D;

         // Loop over first transform
         for(it3D = range3D.first; it3D != range3D.second; ++it3D)
         {
            // Compute third transform
            ForwardConfigurator::integrateND(*it3D, coord);

            range2D = it3D->edgeRange();
            for(it2D = range2D.first; it2D != range2D.second; ++it2D)
            {
               // Compute second transform
               ForwardConfigurator::integrate2D(*it2D, coord);
            }
         }
      } else if(coord.ss().dimension() == 2)
      {
         // Loop over first transform
         for(it3D = range3D.first; it3D != range3D.second; ++it3D)
         {
            // Compute third transform
            ForwardConfigurator::integrateND(*it3D, coord);
         }
      } else
      {
         throw std::logic_error("Configurator cannot be used with less than 2 dimensions");
      }
   }

   template <typename TVariable> void ForwardSingle1DConfigurator::secondStep(const TransformTree&, TVariable&, TransformCoordinatorType&)
   {
      // No need for a second step
   }

   template <typename TVariable> void ForwardSingle1DConfigurator::lastStep(const TransformTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      Profiler::RegionFixture<1> fix("FwdLastStep");

      // Iterators for the transforms
      TransformTreeEdge::EdgeType_citerator it1D;
      TransformTreeEdge::EdgeType_citerator itPhys;

      // Ranges for the vector of edges for the transforms
      TransformTreeEdge::EdgeType_crange range1D;
      TransformTreeEdge::EdgeType_crange rangePhys = tree.root().edgeRange();

      if(coord.ss().dimension() == 3)
      {
         // Iterators for the second transforms
         TransformTreeEdge::EdgeType_citerator it2D;

         // Ranges for the vector of edges for the second transforms
         TransformTreeEdge::EdgeType_crange range2D;

         // Loop over first transform
         for(itPhys = rangePhys.first; itPhys != rangePhys.second; ++itPhys)
         {
            range2D = itPhys->edgeRange();
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
      } else if(coord.ss().dimension() == 2)
      {
         // Loop over first transform
         for(itPhys = rangePhys.first; itPhys != rangePhys.second; ++itPhys)
         {
            range1D = itPhys->edgeRange();
            for(it1D = range1D.first; it1D != range1D.second; ++it1D)
            {
               // Compute third transform
               ForwardConfigurator::integrate1D(*it1D, coord);

               // Update equation
               ForwardConfigurator::updateEquation(*it1D, rVariable, coord);
            }
         }
      } else
      {
         throw std::logic_error("Configurator cannot be used with less than 2 dimensions");
      }
   }

}
}

#endif // QUICC_TRANSFORM_FORWARDSINGLE1DCONFIGURATOR_HPP
