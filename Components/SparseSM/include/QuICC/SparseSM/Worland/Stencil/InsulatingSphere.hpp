/**
 * @file InsulatingSphere.hpp
 * @brief Implementation of the full sphere Worland insulating sphere boundary condition stencil
 */

#ifndef QUICC_SPARSESM_WORLAND_STENCIL_INSULATINGSPHERE_HPP
#define QUICC_SPARSESM_WORLAND_STENCIL_INSULATINGSPHERE_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/IWorlandOperator.hpp"
#include "QuICC/SparseSM/Worland/Stencil/InsulatingSphereDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

   /**
    * @brief Implementation of the full sphere Worland insulating sphere boundary condition stencil
    */
   class InsulatingSphere: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param rows    Number of row
          * @param cols    Number of cols
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param l       Harmonic degree l
          */
         InsulatingSphere(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l);

         /**
          * @brief Destructor
          */
         virtual ~InsulatingSphere() = default;

      protected:

      private:
         /**
          * @brief Build triplet representation of matrix
          *
          * @param list List of triplets (row, col, value)
          */
         void buildTriplets(TripletList_t& list) const final;

         /**
          * @brief Implementation of the diagonals
          */
         std::shared_ptr<InsulatingSphereDiags> mpImpl;
   };

} // Stencil
} // Worland
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_WORLAND_STENCIL_INSULATINGSPHERE_HPP
