/**
 * @file I2Y2SphLapl.hpp
 * @brief Implementation of the I^2 Y^2 spherical laplacian sparse operator, with y = ax + b
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I2Y2SPHLAPL_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I2Y2SPHLAPL_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/ISphericalOperator.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   /**
    * @brief Implementation of the I^2 Y^2 spherical laplacian sparse operator, with y = ax + b
    */
   class I2Y2SphLapl: public ISphericalOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param rows    Number of rows
          * @param cols    Number of cols
          * @param lower   Lower bound
          * @param upper   Upper bound
          * @param l       Harmonic degree l
          */
         I2Y2SphLapl(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const Scalar_t l);

         /**
          * @brief Destructor
          */
         virtual ~I2Y2SphLapl() = default;

      protected:

      private:
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
          * @brief Build triplet representation of matrix
          *
          * @param[out] list containing triplets
          */
         virtual void buildTriplets(TripletList_t& list) const;
   };

} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I2Y2SPHLAPL_HPP
