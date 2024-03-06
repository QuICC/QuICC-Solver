/**
 * @file Id.hpp
 * @brief Implementation of the full sphere Bessel (restricted) identity sparse
 * operator
 */

#ifndef QUICC_SPARSESM_BESSEL_ID_HPP
#define QUICC_SPARSESM_BESSEL_ID_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/IdDiags.hpp"
#include "QuICC/SparseSM/IBesselOperator.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

/**
 * @brief Implementation of the full sphere Bessel (restricted) identity sparse
 * operator
 */
class Id : public IBesselOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of row
    * @param cols    Number of cols
    * @param type    Type of Bessel basis
    * @param l       Harmonic degree l
    * @param s Shift of main diagonal
    */
   Id(const int rows, const int cols, const BesselKind type, const int l,
      const int s = 0);

   /**
    * @brief Destructor
    */
   virtual ~Id() = default;

protected:
private:
   /**
    * @brief Build triplet representation of matrix
    *
    * @param list List of triplets (row, col, value)
    */
   virtual void buildTriplets(TripletList_t& list) const override;

   /**
    * @brief Implementation of the diagonals
    */
   std::shared_ptr<IdDiags> mpImpl;

   /**
    * @brief Shift of main diagonal
    */
   int mShift;
};

} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_BESSEL_ID_HPP
