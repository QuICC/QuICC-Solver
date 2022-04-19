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
#include "QuICC/Framework/Selector/ScalarField.hpp"
#include "QuICC/TransformConfigurators/TransformTree.hpp"
#include "QuICC/TransformConfigurators/ForwardConfigurator.hpp"

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
         ForwardSerialConfigurator() {};

         /**
          * @brief Empty destructor
          */
         virtual ~ForwardSerialConfigurator() {};

      private:
   };

   inline void ForwardSerialConfigurator::setup1DCommunication(const int, TransformCoordinatorType&)
   {
   }

   inline void ForwardSerialConfigurator::setup2DCommunication(const int, TransformCoordinatorType&)
   {
   }

   inline void ForwardSerialConfigurator::initiate1DCommunication(TransformCoordinatorType&)
   {
   }

   inline void ForwardSerialConfigurator::initiate2DCommunication(TransformCoordinatorType&)
   {
   }

   template <typename TVariable> void ForwardSerialConfigurator::firstStep(const TransformTree& tree, TVariable& rVariable, Physical::Kernel::SharedIPhysicalKernel spKernel, TransformCoordinatorType& coord)
   {
      ProfilerMacro_start(Debug::Profiler::FWDTRANSFORM);

      // Iterators for the transforms
      TransformTreeEdge::EdgeType_citerator itSpec;
      TransformTreeEdge::EdgeType_citerator itPhys;

      // Ranges for the vector of edges for the three transforms
      TransformTreeEdge::EdgeType_crange rangeSpec;
      TransformTreeEdge::EdgeType_crange rangePhys = tree.root().edgeRange();

      ProfilerMacro_stop(Debug::Profiler::FWDTRANSFORM);

      // Compute the physical space kernel
      ForwardConfigurator::nonlinearTerm(tree, spKernel, coord);

      ProfilerMacro_start(Debug::Profiler::FWDTRANSFORM);

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
         // Loop over physical transform
         for(itPhys = rangePhys.first; itPhys != rangePhys.second; ++itPhys)
         {
            // Compute transform
            ForwardConfigurator::integrate1ND(*itPhys, coord);

            // Update equation
            ForwardConfigurator::updateEquation(*itPhys, rVariable, coord);
         }
      } else
      {
         throw std::logic_error("Transform with more than 3 dimensions are not implemented");
      }

      ProfilerMacro_stop(Debug::Profiler::FWDTRANSFORM);
   }

   template <typename TVariable> void ForwardSerialConfigurator::secondStep(const TransformTree&, TVariable&, TransformCoordinatorType&)
   {
      // No need for a second step
   }

   template <typename TVariable> void ForwardSerialConfigurator::lastStep(const TransformTree&, TVariable&, TransformCoordinatorType&)
   {
      // No need for a last step
   }

}
}

#endif // QUICC_TRANSFORM_FORWARDSERIALCONFIGURATOR_HPP
