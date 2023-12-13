/**
 * @file IWorlandRadialPower.hpp
 * @brief Interface for a Worland based power operator on radial grid
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_POWER_KOKKOS_IWORLANDRADIALPOWER_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_POWER_KOKKOS_IWORLANDRADIALPOWER_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/IWorlandRadialPower.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/KokkosIWorlandReductor.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

/**
 * @brief Interface for a Worland based POWER
 */
class KokkosIWorlandRadialPower : public KokkosIWorlandReductor
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
   // Apply operator using the entire Matrix in parallel instead  of block by
   // block
   virtual void applyUnitOperator(const OpMatrixL& rOut, const OpMatrixLZL& in,
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
    * @brief Compute power (integral of squared values)
    *
    * @param rOut Output physical values
    * @param in   Input spectral coefficients
    */
   void applyOperators(Matrix& rOut, const MatrixZ& in) const override;

   /**
    * @brief Compute power (integral of squared values)
    *
    * @param rOut Output physical values
    * @param in   Input spectral coefficients
    */
   void applyOperators(MatrixZ& rOut, const MatrixZ& in) const override;
};

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_WORLAND_KOKKOS_IWORLANDRADIALPOWER_HPP
