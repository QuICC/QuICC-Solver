/** 
 * @file I6Diags.hpp
 * @brief Interface to I6 diagonals for full sphere Worland I6 sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_CYLENERGY_I6DIAGS_HPP
#define QUICC_SPARSESM_WORLAND_CYLENERGY_I6DIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Worland/I6Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace CylEnergy {

   /**
    * @brief Implementation of the full sphere Worland I6 sparse operator
    */ 
   class I6Diags: public QuICC::SparseSM::Worland::I6Diags
   {
      public:
         /**
          * @brief Constructor
          *
          * @param alpha   Jacobi alpha
          * @param l       Harmonic degree
          * @param q       Truncation q
          */
         I6Diags(const Scalar_t alpha, const int l, const int q);

         /**
          * @brief Destructor
          */
         virtual ~I6Diags() = default;

         /**
          * @brief 6. subdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d_6(const ACoeff_t& n) const; 

         /**
          * @brief 5. subdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d_5(const ACoeff_t& n) const; 

         /**
          * @brief 4. subdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d_4(const ACoeff_t& n) const; 

         /**
          * @brief 3. subdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d_3(const ACoeff_t& n) const; 

         /**
          * @brief 2. subdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d_2(const ACoeff_t& n) const; 

         /**
          * @brief 1. subdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d_1(const ACoeff_t& n) const; 

         /**
          * @brief Main diagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d0(const ACoeff_t& n) const; 

         /**
          * @brief 1. superdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d1(const ACoeff_t& n) const; 

         /**
          * @brief 2. superdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d2(const ACoeff_t& n) const; 

         /**
          * @brief 3. superdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d3(const ACoeff_t& n) const; 

         /**
          * @brief 4. superdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d4(const ACoeff_t& n) const; 

         /**
          * @brief 5. superdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d5(const ACoeff_t& n) const; 

         /**
          * @brief 6. superdiagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d6(const ACoeff_t& n) const; 
         
      protected:

      private:
   };

}
}
}
}

#endif // QUICC_SPARSESM_WORLAND_CYLENERGY_I6DIAGS_HPP
