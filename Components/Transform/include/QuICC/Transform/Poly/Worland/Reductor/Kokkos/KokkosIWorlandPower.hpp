/**
 * @file IWorlandPower.hpp
 * @brief Interface for a Worland based power operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_POWER_KOKKOS_IWORLANDPOWER_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_POWER_KOKKOS_IWORLANDPOWER_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Kokkos/KokkosIWorlandReductor.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

/**
 * @brief Interface for a Worland based POWER
 */
class KokkosIWorlandPower : public KokkosIWorlandReductor
{
public:
   using OpVectorI = typename KokkosIOperatorTypes::OpVectorI;
   using OpMatrixLZ = typename KokkosIOperatorTypes::OpMatrixLZ;
   using OpMatrixLZL = typename KokkosIOperatorTypes::OpMatrixLZL;
   using OpMatrixL = typename KokkosIOperatorTypes::OpMatrixL;

   /**
    * @brief Constructor
    */
   KokkosIWorlandPower(const int shift);

   /**
    * @brief Destructor
    */
   virtual ~KokkosIWorlandPower();

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
    * @brief Storage for the power integrator
    */
   mutable OpMatrixL vmEOps;

   /**
    * @brief Polynomial shift
    */
   const int mcShift;

   /**
    * @brief Compute power quadrature
    */
   void computePowerQuadrature(Internal::Array& igrid,
      Internal::Array& iweights, const int gSize) const;


   // Apply operator using the entire Matrix in parallel instead  of block by
   // block
   virtual void applyUnitOperator(const OpMatrixL& rOut, const OpMatrixLZL& in,
      const OpVectorI& scan, const int totalOpsCols) const = 0;

   virtual void defaultApplyUnitOperator(const OpMatrixL& rOut,
      const OpMatrixLZL& in, const OpVectorI& scan,
      const int totalOpsCols) const;


private:
   /**
    * @brief Initialise the operators
    */
   virtual void initOperators(const Internal::Array& igrid,
      const Internal::Array& iweights) const override;

   /**
    * @brief Make operator
    */
   virtual void makeOperator(Matrix& op, Matrix& eop,
      const Internal::Array& igrid, const Internal::Array& iweights,
      const int i) const = 0;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_KOKKOS_IWORLANDPOWER_HPP
