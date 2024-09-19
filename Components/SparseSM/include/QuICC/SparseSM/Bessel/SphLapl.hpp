/**
 * @file SphLapl.hpp
 * @brief Implementation of the full sphere Bessel SphLapl sparse operator
 */

#ifndef QUICC_SPARSESM_BESSEL_SPHLAPL_HPP
#define QUICC_SPARSESM_BESSEL_SPHLAPL_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/SphLaplDiags.hpp"
#include "QuICC/SparseSM/IBesselOperator.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

/**
 * @brief Implementation of the full sphere Bessel SphLapl sparse operator
 */
class SphLapl : public IBesselOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of row
    * @param cols    Number of cols
    * @param type    Type of Besse basis
    * @param l       Harmonic degree l
    */
   SphLapl(const int rows, const int cols, const BesselKind type, const int l);

   /**
    * @brief Destructor
    */
   virtual ~SphLapl() = default;

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
   std::shared_ptr<SphLaplDiags> mpImpl;
};

} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_BESSEL_SPHLAPL_HPP
