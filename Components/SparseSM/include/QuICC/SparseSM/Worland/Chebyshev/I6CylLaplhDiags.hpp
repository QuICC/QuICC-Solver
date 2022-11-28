/**
 * @file I6CylLaplhDiags.hpp
 * @brief Interface to I6CylLaplh diagonals for full sphere Worland I6CylLaplh sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_CHEBYSHEV_I6CylLaplhDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_CHEBYSHEV_I6CylLaplhDIAGS_HPP

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
#include "QuICC/SparseSM/Worland/I6CylLaplhDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   /**
    * @brief Implementation of the full sphere Worland I6CylLaplh sparse operator
    */
   class I6CylLaplhDiags: public QuICC::SparseSM::Worland::I6CylLaplhDiags
   {
      public:
         /**
          * @brief Constructor
          */
         I6CylLaplhDiags(const Scalar_t alpha, const int l);

         /**
          * @brief Destructor
          */
         virtual ~I6CylLaplhDiags() = default;

         /**
          * @brief 5. subdiagonal
          */
         ACoeff_t d_5(const ACoeff_t& n) const final;

         /**
          * @brief 4. subdiagonal
          */
         ACoeff_t d_4(const ACoeff_t& n) const final;

         /**
          * @brief 3. subdiagonal
          */
         ACoeff_t d_3(const ACoeff_t& n) const final;

         /**
          * @brief 2. subdiagonal
          */
         ACoeff_t d_2(const ACoeff_t& n) const final;

         /**
          * @brief 1. subdiagonal
          */
         ACoeff_t d_1(const ACoeff_t& n) const final;

         /**
          * @brief Main diagonal
          */
         ACoeff_t d0(const ACoeff_t& n) const final;

         /**
          * @brief 1. superdiagonal
          */
         ACoeff_t d1(const ACoeff_t& n) const final;

         /**
          * @brief 2. superdiagonal
          */
         ACoeff_t d2(const ACoeff_t& n) const final;

         /**
          * @brief 3. superdiagonal
          */
         ACoeff_t d3(const ACoeff_t& n) const final;

         /**
          * @brief 4. superdiagonal
          */
         ACoeff_t d4(const ACoeff_t& n) const final;

         /**
          * @brief 5. superdiagonal
          */
         ACoeff_t d5(const ACoeff_t& n) const final;

      protected:

      private:
   };

}
}
}
}

#endif // QUICC_SPARSESM_WORLAND_CHEBYSHEV_I6CylLaplhDIAGS_HPP
