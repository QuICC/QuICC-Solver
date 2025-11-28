/**
 * @file IBesselOperator.hpp
 * @brief Implementation of the generic interface to the full sphere spherical
 * Bessel sparse operator
 */

#ifndef QUICC_SPARSESM_IBESSELOPERATOR_HPP
#define QUICC_SPARSESM_IBESSELOPERATOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Bessel/BesselKind.hpp"
#include "QuICC/SparseSM/ISparseSMOperator.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

/**
 * @brief Implementation of the generic interface to the full sphere spherical
 * Bessel sparse operator
 */
class IBesselOperator : public ISparseSMOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of columns
    * @param type    Type of Bessel basis
    */
   IBesselOperator(const int rows, const int cols,
      const Bessel::BesselKind type);

   /**
    * @brief Destructor
    */
   virtual ~IBesselOperator() = default;

protected:
   /**
    * @brief Type of Bessel implementation
    */
   Bessel::BesselKind type() const;

private:
   /**
    * @Brief Type of Bessel basis
    */
   Bessel::BesselKind mType;
};

} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_IBESSELOPERATOR_HPP
