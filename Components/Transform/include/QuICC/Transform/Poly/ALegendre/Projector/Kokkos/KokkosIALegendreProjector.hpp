/**
 * @file IALegendreProjector.hpp
 * @brief Interface for a associated Legendre based Projector
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_IALEGENDREPROJECTOR_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_IALEGENDREPROJECTOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/ALegendre/Projector/IALegendreProjector.hpp"
#include "QuICC/Transform/Poly/ALegendre/KokkosIALegendreOperatorTypes.hpp"
#ifdef QUICC_USE_KOKKOS_CUDA
#include "QuICC/Transform/Poly/ALegendre/CudaIALegendreOperatorTypes.hpp"
#include "QuICC/Transform/Poly/ALegendre/CudaIALegendreOperatorGemmUtils.hpp"
#endif
namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Projector {

   /**
    * @brief Interface for a associated Legendre based Projector
    */
   class KokkosIALegendreProjector: public IALegendreProjector
   {
      public:
         using DataType = typename KokkosIALegendreOperatorTypes::DataType;

         using OpVectorI = typename KokkosIALegendreOperatorTypes::OpVectorI;
         using OpMatrixI = typename KokkosIALegendreOperatorTypes::OpMatrixI;
         using OpMatrixLZ = typename KokkosIALegendreOperatorTypes::OpMatrixLZ;
         using OpMatrixL = typename KokkosIALegendreOperatorTypes::OpMatrixL;

         /**
          * @brief Constructor
          */
         KokkosIALegendreProjector() = default;

         /**
          * @brief Destructor
          */
         virtual ~KokkosIALegendreProjector() = default;

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
         mutable OpArray  mGrid;

         /**
          * @brief Storage for the quadrature weights
          */
         mutable OpArray  mWeights;

         //Apply operator using the entire OpMatrix in parallel instead  of block by block

         virtual void applyUnitOperator(const OpMatrixLZ &rOut,
            const OpMatrixLZ &in, const OpVectorI &scan,
            const int totalOpsCols) const = 0;

      private:
         /**
          * @brief Initialise the operators
          */
         virtual void initOperators(const OpArray& igrid, const OpArray& iweights) const;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const = 0;

         /**
          * @brief Compute projection
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void applyOperators(OpMatrixZ &rOut, const OpMatrixZ &in) const;

         /**
          * @brief Apply ith operator
          * this is not used in the Kokkos/Cuda implementations
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

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_PROJECTOR_KOKKOS_IALEGENDREPROJECTOR_HPP
