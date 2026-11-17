/**
 * @file HighamExponential.hpp
 * @brief Exponential of matrix algorithm from Higham 2005
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_HIGHAMEXPONENTIAL_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_HIGHAMEXPONENTIAL_HPP

// System includes
//
#include <array>

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Timestep {

/// @brief This namespace contains exponential timestepping schemes
namespace Exponential {

/**
 * @brief Algorith for exponential of matrix from Higham 2005
 */
class HighamExponential
{
public:
   /**
    * @brief ctor
    */
   HighamExponential();

   /**
    * @brief dtor
    */
   virtual ~HighamExponential() = default;

   /**
    * @brief Compute exponential of matrix
    */
   Matrix compute(const Matrix& mat) const;

private:
   /**
    * @brief Initialize Pade coefficients
    */
   void initPade();

   /**
    * @brief Preprocess A to reduce norm
    *
    * @param scale Diagonal scaling
    * @param ilo   Startindex
    * @param ihi   End index
    * @param matA  Matrix A
    */
   void preprocess(std::vector<double>& scale, int& ilo, int& ihi,
      Matrix& matA) const;

   /**
    * @brief Undo preprocessing
    *
    * @param matA Matrix A
    * @param mu   Diagonal scaling
    * @param matD Balancing matrix
    */
   void postprocess(Matrix& matA, std::vector<double>& scale, int ilo,
      int ihi) const;

   /**
    * @brief Compute matrix norm
    */
   double norm(const Matrix& matA) const;

   /**
    * @brief Theta_m thresholds
    */
   std::array<double, 5> mTheta;

   /**
    * @brief Pade approximant of order 3
    */
   std::array<double, 4> mPade3;

   /**
    * @brief Pade approximant of order 5
    */
   std::array<double, 6> mPade5;

   /**
    * @brief Pade approximant of order 7
    */
   std::array<double, 8> mPade7;

   /**
    * @brief Pade approximant of order 9
    */
   std::array<double, 10> mPade9;

   /**
    * @brief Pade approximant of order 13
    */
   std::array<double, 14> mPade13;
};

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_HIGHAMEXPONENTIAL_HPP
