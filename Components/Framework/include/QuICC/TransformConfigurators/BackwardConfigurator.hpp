/**
 * @file BackwardConfigurator.hpp
 * @brief This class defines the base operations for a backward transform in xD space
 */

#ifndef QUICC_TRANSFORM_BACKWARDCONFIGURATOR_HPP
#define QUICC_TRANSFORM_BACKWARDCONFIGURATOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/TypeSelectors/TransformCommSelector.hpp"
#include "QuICC/TransformConfigurators/TransformTree.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the base operations for a backward transform in xD space
    */
   class BackwardConfigurator
   {
      public:

      protected:
         /**
          * @brief Prepare computation of projection for a scalar variable
          *
          * @param rScalar Scalar variable
          * @param coord   Transform coordinator
          */
         static void prepareSpectral(const TransformTree& tree, Framework::Selector::VariantSharedScalarVariable& rScalar, TransformCoordinatorType& coord);

         /**
          * @brief Prepare computation of projection for a vector variable
          *
          * @param rVector Vector variable
          * @param coord   Transform coordinator
          */
         static void prepareSpectral(const TransformTree& tree, Framework::Selector::VariantSharedVectorVariable& rVector, TransformCoordinatorType& coord);

         /**
          * @brief Prepare computation of physical projection for a scalar variable
          *
          * @param rScalar Scalar variable
          * @param coord   Transform coordinator
          */
         static void preparePhysical(const TransformTree& tree, const TransformTreeEdge& edge, Framework::Selector::VariantSharedScalarVariable& rScalar, TransformCoordinatorType& coord);

         /**
          * @brief Prepare computation of physical projection for a vector variable
          *
          * @param rVector Vector variable
          * @param coord   Transform coordinator
          */
         static void preparePhysical(const TransformTree& tree, const TransformTreeEdge& edge, Framework::Selector::VariantSharedVectorVariable& rVector, TransformCoordinatorType& coord);

         /**
          * @brief Compute the projection transform of the first dimension
          *
          * @param edge    Transform tree edge
          * @param coord   Transform coordinator
          */
         static void project1D(const TransformTreeEdge& edge, TransformCoordinatorType& coord);

         /**
          * @brief Compute the projection transform of an intermediate dimension
          *
          * @param edge    Transform tree edge
          * @param coord   Transform coordinator
          */
         static void project2D(const TransformTreeEdge& edge, TransformCoordinatorType& coord);


         /**
          * @brief Compute the projection transform of the last dimension
          *
          * @param edge    Transform tree edge
          * @param coord   Transform coordinator
          */
         static void projectND(const TransformTreeEdge& edge, TransformCoordinatorType& coord);

         /**
          * @brief Empty constructor
          */
         BackwardConfigurator() {};

         /**
          * @brief Empty destructor
          */
         virtual ~BackwardConfigurator() {};

      private:
         /**
          * @brief Generic implementation of projection (backward transform)
          *
          * @param coord   Transform coordinator
          */
         static void genericProjection(const TransformTreeEdge& edge, TransformCoordinatorType& coord, const Dimensions::Transform::Id, const bool processOutput, const std::string profRegion);
   };

}
}

#endif // QUICC_TRANSFORM_BACKWARDCONFIGURATOR_HPP
