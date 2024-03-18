/**
 * @file CoriolisQm.hpp
 * @brief Implementation of the full sphere Bessel Coriolis cross term acting on l-1
 */

#ifndef QUICC_DENSESM_BESSEL_CORIOLISQM_HPP
#define QUICC_DENSESM_BESSEL_CORIOLISQM_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "DenseSM/IEmbeddedSMOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Bessel {

   /**
    * @brief Implementation of the full sphere Bessel Coriolis cross term acting on l-1
    */
   class CoriolisQm: public IEmbeddedSMOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @paramd outDNu Bessel dNu of output (test function)
          * @paramd inDNu  Bessel dNu of input (basis function)
          * @param rows    Number of row
          * @param cols    Number of cols
          * @param l       Harmonic degree l
          */
         CoriolisQm(const Internal::MHDFloat outDNu, const Internal::MHDFloat inDNu, const int rows, const int cols, const int l);

         /**
          * @brief Destructor
          */
         virtual ~CoriolisQm() = default;

      protected:
         /**
          * @brief Implementation of build dense matrix operator
          *
          * @param mat operator
          * @param rows rows of matrix
          * @param cols cols of matrix
          */
         void buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const final;

         /**
          * @brief Bessel parameter nu = l + dnu for output (test function)
          */
         Internal::MHDFloat mOutDNu;

         /**
          * @brief Bessel parameter nu = l + dnu for input (basis function)
          */
         Internal::MHDFloat mInDNu;

         /**
          * @brief Spherical harmonic degree
          */
         int mL;

      private:
         /**
          * @brief Build generic operator
          *
          * @param rows    Number of row
          * @param cols    Number of cols
          */
         void buildGenericOp(Internal::Matrix& mat, const int rows, const int cols) const;
   };

} // namespace Bessel
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_BESSEL_CORIOLISQM_HPP
