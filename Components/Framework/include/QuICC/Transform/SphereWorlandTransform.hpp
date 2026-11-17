/**
 * @file SphereWorlandTransform.hpp
 * @brief Implementation of the Worland transform in a sphere
 */

#ifndef QUICC_TRANSFORM_SPHEREWORLANDTRANSFORM_HPP
#define QUICC_TRANSFORM_SPHEREWORLANDTRANSFORM_HPP

// System includes
//

// Project includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Transform/ITransform.hpp"
#include "QuICC/Transform/Poly/Worland/Transform.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Implementation of the Worland transform in a sphere
    */
   class SphereWorlandTransform: public ITransform
   {
      public:
         /// Typedef for setup type
         typedef Poly::Worland::Transform::SetupType SetupType;

         /// Typedef for setup type
         typedef Poly::Worland::Transform::SharedSetupType SharedSetupType;

         /**
          * @brief Constructor
          */
         SphereWorlandTransform() = default;

         /**
          * @brief Destructor
          */
         virtual ~SphereWorlandTransform() = default;

         /**
          * @brief set list of required options
          */
         virtual void requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const override;

         /**
          * @brief Set the required options
          */
         virtual void setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId) override;

         /**
          * @brief Get the physical grid
          */
         virtual Array meshGrid() const override;

         /**
          * @brief Initialise the polynomial transform (matrices, weights, grid, etc)
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedSetupType spSetup);

         /**
          * @brief Compute quadrature integration
          *
          * @param rOut       Output spectral coefficients
          * @param in         Input physical values
          * @param id         Integrator to use
          */
         virtual void forward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id) override;

         /**
          * @brief Compute polynomial projection
          *
          * @param rOut       Output physical values
          * @param in         Input spectral coefficients
          * @param id         Projector to use
          */
         virtual void backward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id) override;

         /**
          * @brief Compute reduction operation
          *
          * @param spectrum   Output energy spectrum
          * @param in         Input spectral coefficients
          * @param id         Energy reductor to use
          */
         virtual void reduce(Eigen::Ref<Matrix> spectrum, const Eigen::Ref<const MatrixZ>& in, const std::size_t id) override;

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const override;

         /**
          * @brief Profile the memory requirements
          */
         virtual void profileStorage() const override;

      protected:

      private:
         /**
          * @brief Initialise the operators
          */
         void initOperators();

         /**
          * @brief Transform implementation
          */
         Poly::Worland::Transform mImpl;
   };

} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_SPHEREWORLANDTRANSFORM_HPP
