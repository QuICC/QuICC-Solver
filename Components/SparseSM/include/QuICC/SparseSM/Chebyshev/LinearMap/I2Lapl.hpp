/** 
 * @file I2Lapl.hpp
 * @brief Implementation of the I^2 plane layer cartesian laplaciana operator
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I2LAPL_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I2LAPL_HPP

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
    * @brief Implementation of the I^2 plane layer cartesian laplacian operator
    */ 
   class I2Lapl: public IPlaneOperator
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
         I2Lapl(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const Scalar_t k1, const Scalar_t k2);

         /**
          * @brief Destructor
          */
         virtual ~I2Lapl() = default;
         
      protected:

      private:
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

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I2LAPL_HPP
