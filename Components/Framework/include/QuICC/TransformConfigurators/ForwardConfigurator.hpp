/**
 * @file ForwardConfigurator.hpp
 * @brief This class defines the base operations for a forward transform in xD space
 */

#ifndef QUICC_TRANSFORM_FORWARDCONFIGURATOR_HPP
#define QUICC_TRANSFORM_FORWARDCONFIGURATOR_HPP

// Configuration includes
//
#include "QuICC/Debug/Profiler/ProfilerMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/PhysicalKernels/IPhysicalKernel.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/TypeSelectors/TransformCommSelector.hpp"
#include "QuICC/TransformConfigurators/TransformTree.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the base operations for a forward transform in xD space
    */
   class ForwardConfigurator
   {
      public:

      protected:
         /**
          * @brief Compute nonlinear interaction on a scalar or vector field
          *
          * @param spKernel Physical kernel providing the nonlinear computation
          * @param coord    Transform coordinator
          */
         static void nonlinearTerm(const TransformTree& tree, Physical::Kernel::SharedIPhysicalKernel spKernel, TransformCoordinatorType& coord);

         /**
          * @brief Compute the integration transform of the first dimension
          *
          * @param coord   Transform coordinator
          */
         static void integrate1D(const TransformTreeEdge& edge, TransformCoordinatorType& coord);

         /**
          * @brief Compute the integration transform of the second dimension
          *
          * @param coord   Transform coordinator
          */
         static void integrate2D(const TransformTreeEdge& edge, TransformCoordinatorType& coord);

         /**
          * @brief Compute the integration transform of the last dimension
          *
          * @param coord   Transform coordinator
          */
         static void integrateND(const TransformTreeEdge& edge, TransformCoordinatorType& coord);

         /**
          * @brief Compute the integration transform for the single dimension case
          *
          * @param coord   Transform coordinator
          */
         static void integrate1ND(const TransformTreeEdge& edge, TransformCoordinatorType& coord);

         /**
          * @brief Update scalar variable from dealiased data
          *
          * @param rScalar Scalar variable
          * @param coord   Transform coordinator
          */
         static void updateEquation(const TransformTreeEdge& edge, Framework::Selector::VariantSharedScalarVariable& rScalar, TransformCoordinatorType& coord);

         /**
          * @brief Update vector variable from dealiased data
          *
          * @param rVector Vector variable
          * @param coord   Transform coordinator
          */
         static void updateEquation(const TransformTreeEdge& edge, Framework::Selector::VariantSharedVectorVariable& rVector, TransformCoordinatorType& coord);

         /**
          * @brief Empty constructor
          */
         ForwardConfigurator() {};

         /**
          * @brief Empty destructor
          */
         virtual ~ForwardConfigurator() {};

      private:
   };

}
}

#endif // QUICC_TRANSFORM_FORWARDCONFIGURATOR_HPP
