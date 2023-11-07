/**
 * @file IWorlandIntegrator.hpp
 * @brief Interface for a Worland based integrator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_KOKKOS_IWORLANDINTEGRATOR_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_KOKKOS_IWORLANDINTEGRATOR_HPP


// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Integrator/IWorlandIntegrator.hpp"
#include "QuICC/Transform/Poly/KokkosIOperatorTypes.hpp"
#include "QuICC/Transform/Poly/KokkosIOperatorGemmUtils.hpp"
#ifdef QUICC_USE_KOKKOS_CUDA
#include "QuICC/Transform/Poly/CudaIOperatorTypes.hpp"
#include "QuICC/Transform/Poly/CudaIOperatorGemmUtils.hpp"
#endif

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

/**
 * @brief Interface for a Worland based integrator
 */
class KokkosIWorlandIntegrator : public IWorlandIntegrator
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

protected:
   /**
    * @brief Default implementation to apply ith operator
    */
   void defaultApplyOperator(OpMatrixR rOut, const int i,
      const OpMatrixCR& in) const;

   /**
    * @brief Storage for the operators
    */
   mutable OpMatrixL vmOps;

   /**
    * @brief Storage for the quadrature grid
    */
   mutable Internal::Array mGrid;

   /**
    * @brief Storage for the quadrature weights
    */
   mutable Internal::Array mWeights;

   // Apply operator using the entire Matrix in parallel instead  of block by
   // block

   virtual void applyUnitOperator(const OpMatrixLZ& rOut, const OpMatrixLZL& in,
      const OpVectorI& scan, const int totalOpsCols) const = 0;

private:
   /**
    * @brief Initialise the operators
    */
   virtual void initOperators(const Internal::Array& igrid,
      const Internal::Array& iweights) const override;

   /**
    * @brief Make operator
    */
   virtual void makeOperator(Matrix& op, const Internal::Array& igrid,
      const Internal::Array& iweights, const int i) const = 0;

   /**
    * @brief Compute polynomial integration
    *
    * @param rOut Output spectral coefficients
    * @param in   Input physical values
    */
   void applyOperators(MatrixZ& rPhysVal,
      const MatrixZ& specVal) const override;

   /**
    * @brief Apply ith operator
    */
   void applyOperator(OpMatrixR rOut, const int i,
      const OpMatrixCR& in) const final;
};

} // namespace Integrator
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_WORLAND_IWORLANDINTEGRATOR_HPP
