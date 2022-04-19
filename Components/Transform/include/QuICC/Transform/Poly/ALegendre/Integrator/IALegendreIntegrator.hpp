/** 
 * @file IALegendreIntegrator.hpp
 * @brief Interface for a associated Legendre based integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_IALEGENDREINTEGRATOR_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_IALEGENDREINTEGRATOR_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Poly/ALegendre/IALegendreOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   /**
    * @brief Interface for a associated Legendre based integrator
    */ 
   class IALegendreIntegrator: public IALegendreOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IALegendreIntegrator();

         /**
          * @brief Destructor
          */
         virtual ~IALegendreIntegrator();

         /**
          * @brief Rows of output data
          */
         virtual int outRows() const override;

         /**
          * @brief Columns of output data
          */
         virtual int outCols() const override;

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const override;
         
      protected:
         /**
          * @brief Storage for the operators 
          */
         mutable std::vector<Matrix>  mOps;

         /**
          * @brief Storage for the quadrature grid
          */
         mutable internal::Array  mGrid;

         /**
          * @brief Storage for the quadrature weights
          */
         mutable internal::Array  mWeights;

      private:
         /**
          * @brief Initialise the operators
          */
         virtual void initOperators(const internal::Array& igrid, const internal::Array& iweights) const override;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const = 0;

         /**
          * @brief Compute projection
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void applyOperators(MatrixZ& rOut, const MatrixZ& in) const override;

         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const = 0;

         /**
          * @brief Allow specific initialisation for operators
          */
         virtual void initSpecial() const;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_IALEGENDREINTEGRATOR_HPP
