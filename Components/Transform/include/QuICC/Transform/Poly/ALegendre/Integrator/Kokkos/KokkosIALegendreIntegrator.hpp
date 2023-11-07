/**
 * @file IALegendreIntegrator.hpp
 * @brief Interface for a associated Legendre based integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_IALEGENDREINTEGRATOR_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_IALEGENDREINTEGRATOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Integrator/IALegendreIntegrator.hpp"
#include "QuICC/Transform/Poly/KokkosIOperatorTypes.hpp"
#include "QuICC/Transform/Poly/KokkosIOperatorGemmUtils.hpp"
#ifdef QUICC_USE_KOKKOS_CUDA
#include "QuICC/Transform/Poly/CudaIOperatorTypes.hpp"
#include "QuICC/Transform/Poly/CudaIOperatorGemmUtils.hpp"
#endif

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   /**
    * @brief Interface for a associated Legendre based integrator
    */
   class KokkosIALegendreIntegrator: public IALegendreIntegrator
   {
      public:

         using DataType = typename KokkosIOperatorTypes::DataType;

         using OpVectorI = typename KokkosIOperatorTypes::OpVectorI;
         using OpMatrixI = typename KokkosIOperatorTypes::OpMatrixI;
         using OpMatrixLZ = typename KokkosIOperatorTypes::OpMatrixLZ;
         using OpMatrixLZL = typename KokkosIOperatorTypes::OpMatrixLZL;
         using OpMatrixL = typename KokkosIOperatorTypes::OpMatrixL;

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
         /* MHDFloat requiredStorage() const; */

      protected:
         /**
          * @brief Storage for the operators
          */
         mutable OpMatrixL vmOps;

         /**
          * @brief Storage for the quadrature grid
          */
         mutable Internal::Array  mGrid;

         /**
          * @brief Storage for the quadrature weights
          */
         mutable Internal::Array  mWeights;

         //Apply operator using the entire Matrix in parallel instead  of block by block

         virtual void applyUnitOperator(const OpMatrixLZ &rOut,
            const OpMatrixLZL &in, const OpVectorI &scan,
            const int totalOpsCols) const = 0;

      private:
         /**
          * @brief Initialise the operators
          */
         virtual void initOperators(const Internal::Array& igrid, const Internal::Array& iweights) const;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const = 0;

         /**
          * @brief Compute projection
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void applyOperators(MatrixZ &rOut, const MatrixZ &in) const;

         /**
          * @brief Apply ith operator
          * this is not used in Kokkos/Cuda operators
          */
         void applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const final;

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

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_KOKKOS_IALEGENDREINTEGRATOR_HPP
