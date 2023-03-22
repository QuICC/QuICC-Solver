/** 
 * @file D1.hpp
 * @brief Implementation of the boundary derivative stencil
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_STENCIL_D1_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_STENCIL_D1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Chebyshev/ILinearMapOperator.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Stencil {

   /**
    * @brief Implementation of the boundary derivative stencil
    */ 
   class D1: public ILinearMapOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param rows    Number of rows
          * @param cols    Number of columns
          * @param lower   Lower bound
          * @param upper   Upper bound
          */
         D1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper);

         /**
          * @brief Destructor
          */
         virtual ~D1() = default;
         
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
          * @brief Build triplet representation of matrix
          *
          * @param[out] list containing triplets
          */
         virtual void buildTriplets(TripletList_t& list) const;
   };

} // Stencil
} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_STENCIL_D1_HPP
