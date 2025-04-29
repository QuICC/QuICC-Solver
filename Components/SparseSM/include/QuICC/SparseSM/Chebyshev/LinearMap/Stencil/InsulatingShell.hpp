/**
 * @file InsulatingShell.hpp
 * @brief Implementation of the boundary insulating shell stencil
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_STENCIL_INSULATINGSHELL_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_STENCIL_INSULATINGSHELL_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/ISphericalOperator.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Stencil {

   /**
    * @brief Implementation of the boundary insulating shell stencil
    */
   class InsulatingShell: public ISphericalOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param rows    Number of rows
          * @param cols    Number of columns
          * @param lower   Lower bound
          * @param upper   Upper bound
          * @param l       Harmonic degree l
          */
         InsulatingShell(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const Scalar_t l);

         /**
          * @brief Destructor
          */
         virtual ~InsulatingShell() = default;

      protected:

      private:
         /**
          * @brief 2nd subdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d_2(const ACoeff_t& n) const;

         /**
          * @brief 1st subdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d_1(const ACoeff_t& n) const;

         /**
          * @brief diagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d0(const ACoeff_t& n) const;

         /**
          * @brief Build triplet representation of matrix
          *
          * @param[out] list containing triplets
          */
         virtual void buildTriplets(TripletList_t& list) const;
   };

} // Stencil
} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_STENCIL_INSULATINGSHELL_HPP
