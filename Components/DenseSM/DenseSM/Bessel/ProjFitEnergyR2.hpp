/**
 * @file ProjFitEnergyR2.hpp
 * @brief Implementation of the full sphere Bessel projection operator onto best energy fit
 */

#ifndef QUICC_DENSESM_BESSEL_PROJFITENERGYR2_HPP
#define QUICC_DENSESM_BESSEL_PROJFITENERGYR2_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "DenseSM/IMatrixSMOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Bessel {

   /**
    * @brief Implementation of the full sphere Bessel projection operator onto best energy fit
    */
   class ProjFitEnergyR2: public IMatrixSMOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param outRows Number of rows for output
          * @param bcId    ID of boundary condition
          * @param rows    Number of row
          * @param cols    Number of cols
          * @param l       Harmonic degree l
          * @param q       Truncation q (only consider rows - q equations)
          */
         ProjFitEnergyR2(const int outRows, const std::size_t bcId, const int rows, const int cols, const int l);

         /**
          * @brief Destructor
          */
         virtual ~ProjFitEnergyR2() = default;

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
          * @brief Output truncation
          */
         int mOutRows;

         /**
          * @brief Spherical harmonic degree
          */
         int mL;

         /**
          * @brief Boundary condition
          */
         std::size_t mBcId;

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

#endif // QUICC_DENSESM_BESSEL_PROJFITENERGYR2_HPP
