/**
 * @file IProjCrossOperator.hpp
 * @brief Implementation of the generic projection of cross product A ^ B
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_IPROJCROSSOPERATOR_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_IPROJCROSSOPERATOR_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/ILinearMapOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

/**
 * @brief Implementation of the generic projection of cross product A ^ B
 */
class IProjCrossOperator : public ILinearMapOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of cols
    * @param q       Order of quasi-inverse
    * @param p       Power of radial prefactor
    * @param lOut    Output harmonic degree
    * @param mOut    Output harmonic degree
    * @param lA      harmonic degree of f
    * @param mA      harmonic degree of f
    * @param lB     Input harmonic degree
    * @param mB     Input harmonic degree
    * @param lower   Lower boundary
    * @param upper   Upper boundary
    */
   IProjCrossOperator(const int rows, const int cols, const int q, const int p, const int lOut,
      const int mOut, const int lA, const int mA, const int lB, const int mB,
      const Scalar_t lower, const Scalar_t upper);

   /**
    * @brief Destructor
    */
   virtual ~IProjCrossOperator() = default;

   /**
    * @brief Operator is imaginary?
    */
   bool isImaginary() const;

   /**
    * @brief Operator is exactly zero?
    */
   bool isZero() const;

protected:
   /**
    * @brief Compute Gaunt's integral
    */
   MHDFloat gaunt(const int lA, const int mA, const int lB, const int mB,
      const int lG, const int mG) const;

   /**
    * @brief Compute Elsasser's integral
    */
   MHDFloat elsasser(const int lA, const int mA, const int lB, const int mB,
      const int lG, const int mG) const;

   /**
    * @brief Apply quasi-inverse
    */
   void applyQI(Internal::Matrix& mat) const;

   /**
    * @brief Operator is exactly zero
    */
   bool mIsZero;

   /**
    * @brief Operator is imaginary
    */
   bool mIsImaginary;

   /**
    * @brief Order of quasi-inverse
    */
   int mQ;

   /**
    * @brief Power of radial prefactor
    */
   const int mP;

   /**
    * @brief Harmonic degree of output
    */
   int mLout;

   /**
    * @brief Harmonic order of output
    */
   int mMout;

   /**
    * @brief Harmonic degree of A
    */
   int mLa;

   /**
    * @brief Harmonic order of A
    */
   int mMa;

   /**
    * @brief Harmonic degree of B
    */
   int mLb;

   /**
    * @brief Harmonic order of B
    */
   int mMb;

private:
};

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_IPROJCROSSOPERATOR_HPP
