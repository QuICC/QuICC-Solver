/** 
 * @file I4Diags.hpp
 * @brief Interface to I4 diagonals for full sphere Worland I4 sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_LEGENDRE_I4DIAGS_HPP
#define QUICC_SPARSESM_WORLAND_LEGENDRE_I4DIAGS_HPP

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
#include "QuICC/SparseSM/Worland/I4Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Legendre {

   /**
    * @brief Implementation of the full sphere Worland I4 sparse operator
    */ 
   class I4Diags: public QuICC::SparseSM::Worland::I4Diags
   {
      public:
         /**
          * @brief Constructor
          */
         I4Diags(const Scalar_t alpha, const int l);

         /**
          * @brief Destructor
          */
         virtual ~I4Diags();

         /**
          * @brief 4. subdiagonal
          */
         virtual ACoeff_t d_4(const ACoeff_t& n) const; 

         /**
          * @brief 3. subdiagonal
          */
         virtual ACoeff_t d_3(const ACoeff_t& n) const; 

         /**
          * @brief 2. subdiagonal
          */
         virtual ACoeff_t d_2(const ACoeff_t& n) const; 

         /**
          * @brief 1. subdiagonal
          */
         virtual ACoeff_t d_1(const ACoeff_t& n) const; 

         /**
          * @brief Main diagonal
          */
         virtual ACoeff_t d0(const ACoeff_t& n) const; 

         /**
          * @brief 1. superdiagonal
          */
         virtual ACoeff_t d1(const ACoeff_t& n) const; 

         /**
          * @brief 2. superdiagonal
          */
         virtual ACoeff_t d2(const ACoeff_t& n) const; 

         /**
          * @brief 3. superdiagonal
          */
         virtual ACoeff_t d3(const ACoeff_t& n) const; 

         /**
          * @brief 4. superdiagonal
          */
         virtual ACoeff_t d4(const ACoeff_t& n) const; 
         
      protected:

      private:
   };

}
}
}
}

#endif // QUICC_SPARSESM_WORLAND_LEGENDRE_I4DIAGS_HPP
