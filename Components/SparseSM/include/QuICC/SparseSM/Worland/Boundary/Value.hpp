/**
 * @file Value.hpp
 * @brief Implementation of the bounary value for Worland polynomials
 */

#ifndef QUICC_SPARSESM_WORLAND_BOUNDARY_VALUE_HPP
#define QUICC_SPARSESM_WORLAND_BOUNDARY_VALUE_HPP

// System includes
//

// Project includes
//
#include "QuICC/Precision.hpp"
#include "QuICC/SparseSM/Worland/Boundary/ICondition.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Boundary {

   /**
    * @brief Implementation of the Worland polynomial
    */
   class Value: public ICondition
   {
      public:
         /**
          * @brief Constructor for specific alpha,beta pair
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param l       Harmonic degree l
          */
         Value(const Scalar_t alpha, const Scalar_t dBeta, const int l);

         /**
          * @brief Destructor
          */
         ~Value() = default;

         /**
          * @brief Compute list of boundary values
          *
          * @param maxN       Highest polynomial
          * @param k          Compute value for Jacobi with (alpha+k,beta+k)
          * @param normalized Normalize the values
          */
         ACoeff_t compute(const int maxN, const int k = 0, const bool normalized = true);

      private:
         /**
          * @brief Compute list of boundary values for Chebyshev basis
          *
          * @param k    Compute value for Jacobi with (alpha+k,beta+k)
          * @param maxN Highest polynomial
          */
         ACoeff_t valueChebyshev(const int maxN, const int k = 0);

         /**
          * @brief Compute list of boundary values for Legendre basis
          *
          * @param k    Compute value for Jacobi with (alpha+k,beta+k)
          * @param maxN Highest polynomial
          */
         ACoeff_t valueLegendre(const int maxN, const int k = 0);

         /**
          * @brief Compute list of boundary values for cylindrical energy basis
          *
          * @param k    Compute value for Jacobi with (alpha+k,beta+k)
          * @param maxN Highest polynomial
          */
         ACoeff_t valueCylEnergy(const int maxN, const int k = 0);

         /**
          * @brief Compute list of boundary values for spherical energy basis
          *
          * @param k    Compute value for Jacobi with (alpha+k,beta+k)
          * @param maxN Highest polynomial
          */
         ACoeff_t valueSphEnergy(const int maxN, const int k = 0);

   };

} // Boundary
} // Worland
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_WORLAND_BOUNDARY_VALUE_HPP
