/**
 * @file Operators.cpp
 * @brief Source of the generic utils
 */

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "DenseSM/Chebyshev/LinearMap/Utils/Operators.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Utils {

Matrix computeExpansion(const Matrix& f, const int fN, const Internal::MHDFloat lb, const Internal::MHDFloat ub)
{
   namespace cheb = Transform::Fft::Chebyshev::LinearMap;
   typedef cheb::Integrator::P::SetupType SetupType;

   const auto pId = GridPurpose::SIMULATION;
   int rN = f.rows();
   int cols = f.cols();

   // Setup integrator P
   auto sFFwd = std::make_shared<SetupType>(rN, cols, fN, pId);
   sFFwd->setBounds(static_cast<MHDFloat>(lb),
      static_cast<MHDFloat>(ub));
   sFFwd->lock();
   cheb::Integrator::P TFFwd;
   TFFwd.init(sFFwd);
   
   // Compute spectral expansion
   Matrix sF = Matrix::Zero(rN, cols);
   TFFwd.transform(sF, f);

   return sF;
}

Matrix evaluate(const Matrix& sf, const int fN, const Internal::MHDFloat lb, const Internal::MHDFloat ub)
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
   cheb::Projector::P TFdBwd;
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
