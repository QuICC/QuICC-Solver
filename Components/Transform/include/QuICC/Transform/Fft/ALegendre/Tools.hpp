/**
 * @file Tools.hpp
 * @brief Tools specific to FFT based ALegendre transform
 */

#ifndef QUICC_TRANSFORM_FFT_ALEGENDRE_TOOLS_HPP
#define QUICC_TRANSFORM_FFT_ALEGENDRE_TOOLS_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace ALegendre {

namespace Tools {

/**
* @brief Compute physical grid
*/
void computeGrid(Internal::Array& grid, const int size);

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_TOOLS_HPP
