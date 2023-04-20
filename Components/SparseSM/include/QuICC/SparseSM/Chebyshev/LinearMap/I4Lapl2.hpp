/** 
 * @file I4Lapl2.hpp
 * @brief Implementation of the I^4  plane layer bilaplacian sparse operator, with y = ax + b
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I4LAPL2_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I4LAPL2_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/IPlaneOperator.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   /**
    * @brief Implementation of the I^4 plane layer bilaplacian sparse operator, with y = ax + b
    */ 
   class I4Lapl2: public IPlaneOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param rows    Number of rows
          * @param cols    Number of columns
          * @param lower   Lower bound of domain
          * @param upper   Upper bound of domain
          * @param k1      Fourier mode in first direction
          * @param k2      Fourier mode in second direction
          */
         I4Lapl2(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const Scalar_t k1, const Scalar_t k2);

         /**
          * @brief Destructor
          */
         virtual ~I4Lapl2() = default;
         
      protected:

      private:
         /**
          * @brief 4th subdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d_4(const ACoeff_t& n) const; 

         /**
          * @brief 2nd subdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d_2(const ACoeff_t& n) const; 

         /**
          * @brief diagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d0(const ACoeff_t& n) const; 

         /**
          * @brief 2nd superdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d2(const ACoeff_t& n) const; 

         /**
          * @brief 4th superdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d4(const ACoeff_t& n) const; 

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

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I4LAPL2_HPP
