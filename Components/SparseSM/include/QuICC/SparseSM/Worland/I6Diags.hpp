/** 
 * @file I6Diags.hpp
 * @brief Interface to I6 diagonals for full sphere Worland I6 sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_I6DIAGS_HPP
#define QUICC_SPARSESM_WORLAND_I6DIAGS_HPP

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
#include "QuICC/SparseSM/Worland/IDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   /**
    * @brief Implementation of the full sphere Worland I6 sparse operator
    */ 
   class I6Diags: public IDiags
   {
      public:
         /**
          * @brief Constructor
          */
         I6Diags(const Scalar_t alpha, const Scalar_t dBeta, const int l);

         /**
          * @brief Destructor
          */
         virtual ~I6Diags();

         /**
          * @brief 6. subdiagonal
          */
         virtual ACoeff_t d_6(const ACoeff_t& n) const = 0; 

         /**
          * @brief 5. subdiagonal
          */
         virtual ACoeff_t d_5(const ACoeff_t& n) const = 0; 

         /**
          * @brief 4. subdiagonal
          */
         virtual ACoeff_t d_4(const ACoeff_t& n) const = 0; 

         /**
          * @brief 3. subdiagonal
          */
         virtual ACoeff_t d_3(const ACoeff_t& n) const = 0; 

         /**
          * @brief 2. subdiagonal
          */
         virtual ACoeff_t d_2(const ACoeff_t& n) const = 0; 

         /**
          * @brief 1. subdiagonal
          */
         virtual ACoeff_t d_1(const ACoeff_t& n) const = 0; 

         /**
          * @brief Main diagonal
          */
         virtual ACoeff_t d0(const ACoeff_t& n) const = 0; 

         /**
          * @brief 1. superdiagonal
          */
         virtual ACoeff_t d1(const ACoeff_t& n) const = 0; 

         /**
          * @brief 2. superdiagonal
          */
         virtual ACoeff_t d2(const ACoeff_t& n) const = 0; 

         /**
          * @brief 3. superdiagonal
          */
         virtual ACoeff_t d3(const ACoeff_t& n) const = 0; 

         /**
          * @brief 4. superdiagonal
          */
         virtual ACoeff_t d4(const ACoeff_t& n) const = 0; 

         /**
          * @brief 5. superdiagonal
          */
         virtual ACoeff_t d5(const ACoeff_t& n) const = 0; 

         /**
          * @brief 6. superdiagonal
          */
         virtual ACoeff_t d6(const ACoeff_t& n) const = 0; 
         
      protected:

      private:
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_I6DIAGS_HPP
