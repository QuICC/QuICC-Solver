/**
 * @file ProjFitEnergyPol.hpp
 * @brief Implementation of the full sphere Worland projection operator onto best energy fit for poloidal scalar
 */

#ifndef QUICC_DENSESM_WORLAND_PROJFITENERGYPOL_HPP
#define QUICC_DENSESM_WORLAND_PROJFITENERGYPOL_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "DenseSM/Worland/IEmbeddedOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   /**
    * @brief Implementation of the full sphere Worland projection operator onto best energy fit for poloidal scalar
    */
   class ProjFitEnergyPol: public IEmbeddedOperator
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
         ProjFitEnergyPol(const int outRows, const std::size_t bcId, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q = 0);

         /**
          * @brief Destructor
          */
         virtual ~ProjFitEnergyPol() = default;

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
          * @brief Build Chebyshev type operator
          */
         void buildChebyshevOp(Internal::Matrix& mat, const int rows, const int cols) const;
   };

} // Worland
} // DenseSM
} // QuICC

#endif // QUICC_DENSESM_WORLAND_PROJFITENERGYPOL_HPP
