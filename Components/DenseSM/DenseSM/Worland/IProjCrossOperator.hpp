/**
 * @file IProjCrossOperator.hpp
 * @brief Implementation of the generic projection of cross product A ^ B
 */

#ifndef QUICC_DENSESM_WORLAND_IPROJCROSSOPERATOR_HPP
#define QUICC_DENSESM_WORLAND_IPROJCROSSOPERATOR_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Worland/IWorlandOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

/**
 * @brief Implementation of the generic projection of cross product A ^ B
 */
class IProjCrossOperator : public IWorlandOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of cols
    * @param q       Order of quasi-inverse
    * @param lOut    Output harmonic degree
    * @param mOut    Output harmonic degree
    * @param lA      harmonic degree of f
    * @param mA      harmonic degree of f
    * @param lB      Input harmonic degree
    * @param mB      Input harmonic degree
    * @param alpha   Jacobi alpha parameter
    * @param dBeta   Jacobi dBeta parameter
    */
   IProjCrossOperator(const int rows, const int cols, const int q, const int lOut,
      const int mOut, const int lA, const int mA, const int lB, const int mB,
      const Scalar_t alpha, const Scalar_t dBeta);

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
    * @brief Apply quasi-inverse
    *
    * @param mat     Operator matrix to apply QI on
    * @param l       harmonic degree
    */
   void applyQI(Internal::Matrix& mat, const int l) const;

   /**
    * @brief Set operator band
    */
   void setBand(const int dL, int s);

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

   /**
    * @brief Bandwidth of operator
    */
   std::pair<int,int> mBand;

private:
};

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_IPROJCROSSOPERATOR_HPP
