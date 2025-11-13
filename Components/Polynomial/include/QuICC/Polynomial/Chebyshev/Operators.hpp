/**
 * @file Tools.hpp
 * @brief Tools specific to Chebyshev polynomial implementation
 */

#ifndef QUICC_POLYNOMIAL_CHEBYSHEV_OPERATORS_HPP
#define QUICC_POLYNOMIAL_CHEBYSHEV_OPERATORS_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Polynomial {

namespace Chebyshev {

namespace Operators {

/**
 * @brief normalization factor for the Chebyshev polynomials
 */
static int fT(int n);  

/**
 * @brief factor to half the first element of the sum in (2) of Saibaba, 2021
 */
static MHDFloat fc(int n);    

/**
 * @brief Integral of Chebyshev polynomial, Tn over -1<=x<=1
 */
static MHDFloat IntT(int n);  

/**
 * @brief Integral of Tn Tj over -1<=x<=1
 */
static MHDFloat IntTnTj(int n, int j);  

/**
 * @brief ratio  IntTnTj(int n, int j-2) / IntTnTj(int n, int j)
 */
static MHDFloat ratioIntTnTj2(int n, int j);  


/**
 * @brief ratio of cBar(n,j-2)/cBar(n,j)
 */
static MHDFloat ratioC2(int n, int j);

/**
 * @brief ratio for the recursive relation to calculate tmpSum in Bmat
 */
static MHDFloat ratioS(int p, int k, int n0, int n, int alpha); 

/**
 * @brief Integrate r^p Tn over r
 */
void integrateRpTn(Matrix& iop, const int p,
   const int size, const MHDFloat ro, const MHDFloat ri);

} // namespace Operators
} // namespace Chebyshev
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_CHEBYSHEV_OPERATORS_HPP
