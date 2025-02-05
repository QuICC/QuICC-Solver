/**
 * @file FR3D3.hpp
 * @brief Implementation of the spectral operator f r^3 D3(*)
 * 
 * Modified from FR2D2.cpp
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_FR3D3_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_FR3D3_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/ITripleHarmonicOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

/**
 * @brief Implementation of the spectral operator f r^3 D3(*)
 */
class FR3D3 : public ITripleHarmonicOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of cols
    * @param lOut    Output harmonic degree
    * @param mOut    Output harmonic order
    * @param lF      harmonic degree of f
    * @param mF      harmonic order of f
    * @param lIn     Input harmonic degree
    * @param mIn     Input harmonic order
    * @param lower   Lower boundary
    * @param pF      pointer to a given function
    * @param upper   Upper boundary
    */
   FR3D3(const int rows, 
         const int cols, 
         const int lOut, 
         const int mOut, 
         const int lF, 
         const int mF, 
         const int lIn, 
         const int mIn, 
         std::shared_ptr<RadialTorPolFunction> pF, 
         const Scalar_t lower,
         const Scalar_t upper);

   /**
    * @brief Destructor
    */
   virtual ~FR3D3() = default;

protected:
   /**
    * @brief Implementation of build dense matrix operator
    *
    * @param mat operator
    * @param rows rows of matrix
    * @param cols cols of matrix
    */
   void buildOpImpl(Internal::Matrix& mat, const int rows,
      const int cols) const final;

private:
};

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_FR3D3_HPP
