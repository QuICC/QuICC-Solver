/**
 * @file Tools.hpp
 * @brief Tools for dense full sphere Worland operators
 */

#ifndef QUICC_DENSESM_WORLAND_TOOLS_HPP
#define QUICC_DENSESM_WORLAND_TOOLS_HPP

// Project includes
//
#include "DenseSM/Worland/WorlandKind.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

/**
 * @brief Tools for dense full sphere Worland operators
 */
class Tools
{
public:
   /// Typedef for scalar
   typedef Internal::MHDFloat Scalar_t;

   /**
    * @brief Identify type of Worland basis
    *
    * @param alpha   Jacobi alpha
    * @param dBeta   Jacobi beta = l + dBeta
    */
   static WorlandKind identifyBasis(const Scalar_t alpha, const Scalar_t dBeta);

protected:
private:
   /**
    * @brief Constructor
    */
   Tools() = default;

   /**
    * @brief Destructor
    */
   virtual ~Tools() = default;
};

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_TOOLS_HPP
