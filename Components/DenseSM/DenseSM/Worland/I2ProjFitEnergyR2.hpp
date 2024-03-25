/**
 * @file I2ProjFitEnergyR2.hpp
 * @brief Implementation of the full sphere Worland projection operator onto best energy fit with I2 quas-inverse
 *
 */

#ifndef QUICC_DENSESM_WORLAND_I2PROJFITENERGYR2_HPP
#define QUICC_DENSESM_WORLAND_I2PROJFITENERGYR2_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "DenseSM/Worland/IWorlandOperator.hpp"
#include "DenseSM/Worland/ProjFitEnergyR2.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   /**
    * @brief Implementation of the full sphere Worland projection operator onto best energy fit with I2 quasi-inverse
    */
   class I2ProjFitEnergyR2: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param outRows Number of rows for output
          * @param bcId    ID of boundary condition
          * @param rows    Number of row
          * @param cols    Number of cols
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param l       Harmonic degree l
          * @param q       Truncation q (only consider rows - q equations)
          */
         I2ProjFitEnergyR2(const int outRows, const std::size_t bcId, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q = 0);

         /**
          * @brief Destructor
          */
         virtual ~I2ProjFitEnergyR2() = default;

      protected:
         /**
          * @brief Implementation of build dense matrix operator
          *
          * @param mat operator
          * @param rows rows of matrix
          * @param cols cols of matrix
          */
         void buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const final;

      private:
         /**
          * @brief Build generic operator
          *
          * @param mat     Input/Output matrix
          * @param rows    Number of row
          * @param cols    Number of cols
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          */
         void buildGenericOp(Internal::Matrix& mat, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta) const;

         /**
          * @brief Output truncation
          */
         int mOutRows;

         /**
          * @brief Spherical harmonic degree
          */
         int mL;

         /**
          * @Brief Energy fit projector
          */
         ProjFitEnergyR2 mProj;
   };

} // Worland
} // DenseSM
} // QuICC

#endif // QUICC_DENSESM_WORLAND_I2PROJFITENERGYR2_HPP
