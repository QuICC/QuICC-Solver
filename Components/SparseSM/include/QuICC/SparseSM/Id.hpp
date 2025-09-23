/**
 * @file Id.hpp
 * @brief Implementation of the restricted identity sparse operator
 */

#ifndef QUICC_SPARSESM_ID_HPP
#define QUICC_SPARSESM_ID_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/ISparseSMOperator.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

/**
 * @brief Implementation of restricted identity sparse operator
 */
class Id : public ISparseSMOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of row
    * @param cols    Number of cols
    * @param q       Truncation q (only consider rows - q equations)
    * @param s Shift of main diagonal
    */
   Id(const int rows, const int cols, const int q = 0, const int s = 0);

   /**
    * @brief Destructor
    */
   virtual ~Id() = default;

protected:
private:
   /**
    * @brief Compute diagonal
    */
   ACoeff_t d(const ACoeff_t& n) const;

   /**
    * @brief Build triplet representation of matrix
    *
    * @param list List of triplets (row, col, value)
    */
   virtual void buildTriplets(TripletList_t& list) const override;

   /**
    * @brief Zero rows at top or bottom
    */
   int mQ;

   /**
    * @brief Shift of main diagonal
    */
   int mShift;
};

} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_ID_HPP
