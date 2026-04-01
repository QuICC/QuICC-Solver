/**
 * @file ALegendreTransform.hpp
 * @brief Implementation of the associated Legendre transform
 */

#ifndef QUICC_TRANSFORM_ALEGENDRETRANSFORM_HPP
#define QUICC_TRANSFORM_ALEGENDRETRANSFORM_HPP

// System includes
//
#include <set>
#include <map>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Transform/ITransform.hpp"
#include "QuICC/Transform/Poly/ALegendre/Transform.hpp"
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

namespace QuICC {

namespace Transform {

   /**
    * @brief Implementation of the associated Legendre transform
    */
   class ALegendreTransform: public ITransform
   {
      public:
         /// Typedef for setup type
         typedef Poly::ALegendre::Transform::SetupType SetupType;

         /// Typedef for setup type
         typedef Poly::ALegendre::Transform::SharedSetupType SharedSetupType;

         /**
          * @brief Constructor
          */
         ALegendreTransform() = default;

         /**
          * @brief Destructor
          */
         ~ALegendreTransform() = default;

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
          * @param rSpecVal   Output spectral coefficients
          * @param physVal    Input physical values
          * @param integrator Integrator to use
          */
         void forward(Eigen::Ref<MatrixZ> rSpecVal, const Eigen::Ref<const MatrixZ>& physVal, const std::size_t integrator) override;

         /**
          * @brief Compute polynomial projection
          *
          * @param rPhysVal   Output physical values
          * @param specVal    Input spectral coefficients
          * @param projector  Projector to use
          */
         void backward(Eigen::Ref<MatrixZ> rPhysVal, const Eigen::Ref<const MatrixZ>& specVal, const std::size_t projector) override;

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
         Poly::ALegendre::Transform mImpl;
   };

}
}

#endif // QUICC_TRANSFORM_ALEGENDRETRANSFORM_HPP
