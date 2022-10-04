/** 
 * @file StorageKind.hpp
 * @brief Enum to select storage kind
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_STORAGEKIND_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_STORAGEKIND_HPP

// System includes
//

// Project includes
//

namespace QuICC {

namespace Transform {

/// This namespace provides all the FFT related code.
namespace Fft {


   /**
    * @brief Selector to get temporary storage
    *
    */
   enum class StorageKind
   {
      in,  ///< is used for padded FFT input
      out, ///< is used for FFT output
      mid, ///< is used for physical space intermediate computations
   };

} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif
