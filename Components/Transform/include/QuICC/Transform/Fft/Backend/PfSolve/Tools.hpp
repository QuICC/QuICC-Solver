/**
 * @file Tools.hpp
 * @brief Definition of some useful constants and tools for PFSOLVE
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_PFSOLVE_TOOLS_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_PFSOLVE_TOOLS_HPP

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace PfSolve {

/**
 * @brief Contains some useful constants and tools for PFSOLVE
 */
class Tools {
public:
  /**
   * @brief Optimise the FFT sizes
   *
   * @param size Current size of the FFT
   * @param opt Optimized size extension
   */
  static bool optimizeFft(const int size, int &opt);

  /**
   * @brief Standard dealiasing factor (usually 3/2)
   */
  static const MHDFloat STD_DEALIASING;

  /**
   * @brief The real <-> complex fast Fourrier transform dealiasing factor
   */
  static const MHDFloat MIXED_DEALIASING;

  /**
   * @brief Cosine dealiasing factor (usually 3/2)
   */
  static const MHDFloat COS_DEALIASING;

  /**
   * @brief Maximul extension width to consider for optimization
   */
  static const MHDFloat OPTIMIZATION_WIDTH;

protected:
private:
  /**
   * @brief Empty constructor
   */
  Tools();

  /**
   * @brief Empty Destructor
   */
  ~Tools();
};

} // namespace PfSolve
} // namespace Backend
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_BACKEND_PFSOLVE_TOOLS_HPP
