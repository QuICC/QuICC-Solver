/**
 * @file I2Y1.hpp
 * @brief Implementation of the I^2 Y sparse operator, with y = ax + b
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I2Y1_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I2Y1_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SparseSM/Chebyshev/ILinearMapOperator.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   /**
    * @brief Implementation of the I^2 Y sparse operator, with y = ax + b
    */
   class I2Y1: public ILinearMapOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param rows    Number of rows
          * @param cols    Number of columns
          * @param lower   Lower bound of domain
          * @param upper   Upper bound of domain
          */
         I2Y1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper);

         /**
          * @brief Destructor
          */
         virtual ~I2Y1() = default;

      protected:

      private:
         /**
          * @brief 3rd subdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d_3(const ACoeff_t& n) const;

         /**
          * @brief 2nd subdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d_2(const ACoeff_t& n) const;

         /**
          * @brief 1st subdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d_1(const ACoeff_t& n) const;

         /**
          * @brief diagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d0(const ACoeff_t& n) const;

         /**
          * @brief 1st superdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d1(const ACoeff_t& n) const;

         /**
          * @brief 2nd superdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d2(const ACoeff_t& n) const;

         /**
          * @brief 3rd superdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d3(const ACoeff_t& n) const;

         /**
          * @brief Build triplet representation of matrix
          *
          * @param[out] list containing triplets
          */
         virtual void buildTriplets(TripletList_t& list) const;
   };

}
}
}
}

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I2Y1_HPP
