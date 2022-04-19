/** 
 * @file I4Y3D1Y1.hpp
 * @brief Implementation of the I^4 Y^3 D Y sparse operator, with y = ax + b
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I4Y3D1Y1_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I4Y3D1Y1_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Chebyshev/ILinearMapOperator.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   /**
    * @brief Implementation of the I^4 Y^3 D Y sparse operator, with y = ax + b
    */ 
   class I4Y3D1Y1: public ILinearMapOperator
   {
      public:
         /**
          * @brief Constructor
          */
         I4Y3D1Y1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper);

         /**
          * @brief Destructor
          */
         virtual ~I4Y3D1Y1();
         
      protected:

      private:
         /**
          * @brief 7th subdiagonal
          */
         ACoeff_t d_7(const ACoeff_t& n) const; 

         /**
          * @brief 6th subdiagonal
          */
         ACoeff_t d_6(const ACoeff_t& n) const; 

         /**
          * @brief 5th subdiagonal
          */
         ACoeff_t d_5(const ACoeff_t& n) const; 

         /**
          * @brief 4th subdiagonal
          */
         ACoeff_t d_4(const ACoeff_t& n) const; 

         /**
          * @brief 3rd subdiagonal
          */
         ACoeff_t d_3(const ACoeff_t& n) const; 

         /**
          * @brief 2nd subdiagonal
          */
         ACoeff_t d_2(const ACoeff_t& n) const; 

         /**
          * @brief 1st subdiagonal
          */
         ACoeff_t d_1(const ACoeff_t& n) const; 

         /**
          * @brief diagonal
          */
         ACoeff_t d0(const ACoeff_t& n) const; 

         /**
          * @brief 1st superdiagonal
          */
         ACoeff_t d1(const ACoeff_t& n) const; 

         /**
          * @brief 2nd superdiagonal
          */
         ACoeff_t d2(const ACoeff_t& n) const; 

         /**
          * @brief 3rd superdiagonal
          */
         ACoeff_t d3(const ACoeff_t& n) const; 

         /**
          * @brief 4th superdiagonal
          */
         ACoeff_t d4(const ACoeff_t& n) const; 

         /**
          * @brief 5th superdiagonal
          */
         ACoeff_t d5(const ACoeff_t& n) const; 

         /**
          * @brief 6th superdiagonal
          */
         ACoeff_t d6(const ACoeff_t& n) const; 

         /**
          * @brief 7th superdiagonal
          */
         ACoeff_t d7(const ACoeff_t& n) const; 

         /**
          * @brief Build triplet representation of matrix
          */
         virtual void buildTriplets(TripletList_t& list) const;
   };

}
}
}
}

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_I4Y3D1Y1_HPP
