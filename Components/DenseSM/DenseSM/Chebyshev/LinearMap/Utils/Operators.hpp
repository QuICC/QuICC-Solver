/**
 * @file Operators.hpp
 * @brief Implementation of the generic spectral triple harmonic operator
 */

#ifndef QUICC_DENSESM_CHEBYSHEV_LINEARMAP_UTILS_OPERATORS_HPP
#define QUICC_DENSESM_CHEBYSHEV_LINEARMAP_UTILS_OPERATORS_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Utils {

/**
 * @brief Compute operator for product with spectral expansion
 *
 * @param mat   Matrix to store operator
 * @param rows  Rows of operator
 * @param cols  Cols of operator
 * @param spec  Spectral expansion coefficients
 * @param spec  Size of expansion
 */
template <typename TMat, typename TSpec> void expansionProduct(TMat& mat, const int rows, const int cols, const TSpec& spec, const int nN);

/**
 * @brief Compute spectral expansion of grid function
 */
Matrix computeExpansion(const Matrix& f, const int fN, const Internal::MHDFloat lb, const Internal::MHDFloat ub);

/**
 * @brief Compute values on grid
 */
Matrix evaluate(const Matrix& sf, const int fN, const Internal::MHDFloat lb, const Internal::MHDFloat ub);

/**
 * @brief Compute spectral derivative on grid
 */
template<int N> Matrix evaluateD(const Matrix& sf, const int fN, const Internal::MHDFloat lb, const Internal::MHDFloat ub);

/**
 * @brief Compute spectral operator on grid
 */
template<typename TOp> Matrix evaluateOp(const Matrix& sf, const int fN, const Internal::MHDFloat lb, const Internal::MHDFloat ub);

template <typename TMat, typename TSpec> void expansionProduct(TMat& mat, const int rows, const int cols, const TSpec& spec, const int nN)
{
   assert(mat.rows() >= rows);
   assert(mat.cols() >= cols);
   assert(spec.rows() >= nN);
   assert(spec.cols() >= 1);

   mat.setZero();
   for(int i = 0; i < rows; i++)
   {
      for(int j = 0; j < nN; j++)
      {
         if(i-j >= 0)
         {
            mat(i,i-j) += spec(j);
         }
         else if(j-i < nN)
         {
            mat(i,j-i) += spec(j);
         }
         if(j > 0 && i+j < cols)
         {
            mat(i,i+j) += spec(j);
         }
      }
   }
}

template<int N> Matrix evaluateD(const Matrix& sf, const int fN, const Internal::MHDFloat lb, const Internal::MHDFloat ub)
{
   namespace cheb = Transform::Fft::Chebyshev::LinearMap;
   typedef cheb::Integrator::P::SetupType SetupType;

   const auto pId = GridPurpose::SIMULATION;
   int rN = sf.rows();
   int cols = sf.cols();

   assert(f.rows() == rN);

   // Setup differentiation projector
   auto sFBwd = std::make_shared<SetupType>(rN, cols, fN, pId);
   sFBwd->setBounds(static_cast<MHDFloat>(lb),
      static_cast<MHDFloat>(ub));
   sFBwd->lock();
   cheb::Projector::D<N> TFdBwd;
   TFdBwd.init(sFBwd);

   Matrix df = Matrix::Zero(rN, cols);
   TFdBwd.transform(df, sf);

   return df;
}

template<typename TOp> Matrix evaluateOp(const Matrix& sf, const int fN, const Internal::MHDFloat lb, const Internal::MHDFloat ub)
{
   namespace cheb = Transform::Fft::Chebyshev::LinearMap;
   typedef cheb::Integrator::P::SetupType SetupType;

   const auto pId = GridPurpose::SIMULATION;
   int rN = sf.rows();
   int cols = sf.cols();

   assert(f.rows() == rN);

   // Setup differentiation projector
   auto sFBwd = std::make_shared<SetupType>(rN, cols, fN, pId);
   sFBwd->setBounds(static_cast<MHDFloat>(lb),
      static_cast<MHDFloat>(ub));
   sFBwd->lock();
   TOp TFdBwd;
   TFdBwd.init(sFBwd);

   Matrix f = Matrix::Zero(rN, cols);
   TFdBwd.transform(f, sf);

   return f;
}

} // namespace Utils
} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_CHEBYSHEV_LINEARMAP_UTILS_OPERATORS_HPP
