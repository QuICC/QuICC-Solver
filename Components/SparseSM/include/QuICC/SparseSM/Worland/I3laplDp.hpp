/**
 * @file I3laplDp.hpp
 * @brief Implementation of the full sphere Worland I3 of lapl*D on l+1 coupling
 * sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_I3LAPLDP_HPP
#define QUICC_SPARSESM_WORLAND_I3LAPLDP_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/IWorlandOperator.hpp"
#include "QuICC/SparseSM/Worland/I3laplDpDiags.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

/// Namespace for sparse spectral method operators
namespace SparseSM {

/// Namespace for sparse worland operators
namespace Worland {

/**
 * @brief Implementation of the full sphere Worland I3 of lapl*D on l+1 coupling
 * sparse operator (Stage 4-5 spheroid LHS; never truncated, q>0 unsupported)
 */
class I3laplDp : public IWorlandOperator
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
    * @param q       Truncation q (only consider rows - q equations)
    */
   I3laplDp(const int rows, const int cols, const Scalar_t alpha,
      const Scalar_t dBeta, const int l, const int q = 0);

   /**
    * @brief Destructor
    */
   virtual ~I3laplDp() = default;

protected:
private:
   /**
    * @brief Build triplet representation of matrix
    *
    * @param list List of triplets (row, col, value)
    */
   void buildTriplets(TripletList_t& list) const final;

   /**
    * @brief Build BLAS banded representation of matrix
    *
    * @param bd   Matrix entries in BLAS banded format
    * @param kL   Number of lower diagonals
    * @param kU   Number of upper diagonals
    */
   void buildBanded(Internal::Matrix& bd, unsigned int& kL,
      unsigned int& kU) const final;

   /**
    * @brief Implementation of the diagonals
    */
   std::shared_ptr<I3laplDpDiags> mpImpl;
};

} // namespace Worland
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_WORLAND_I3LAPLDP_HPP
