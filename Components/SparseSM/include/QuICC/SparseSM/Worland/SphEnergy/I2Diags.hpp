/** 
 * @file I2Diags.hpp
 * @brief Interface to I2 diagonals for full sphere Worland I2 sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_SPHENERGY_I2DIAGS_HPP
#define QUICC_SPARSESM_WORLAND_SPHENERGY_I2DIAGS_HPP

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
#include "QuICC/SparseSM/Worland/I2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace SphEnergy {

   /**
    * @brief Implementation of the full sphere Worland I2 sparse operator
    */ 
   class I2Diags: public QuICC::SparseSM::Worland::I2Diags
   {
      public:
         /**
          * @brief Constructor
          */
         I2Diags(const Scalar_t alpha, const int l);

         /**
          * @brief Destructor
          */
         virtual ~I2Diags();

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
         
      protected:

      private:
   };

}
}
}
}

#endif // QUICC_SPARSESM_WORLAND_SPHENERGY_I2DIAGS_HPP
