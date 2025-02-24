/**
 * @file Operators.cpp
 * @brief Source of the generic utils
 */

// System includes
//

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/Utils/Operators.hpp"
#include "Types/Typedefs.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Id.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I3.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4.hpp"

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

   assert(sf.rows() == rN);

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

SparseMatrix selectIq(const int q, const int rows, const int cols, const Internal::MHDFloat lb, const Internal::MHDFloat ub)
{
   if(q < 1)
   {
      throw std::logic_error("Cannot have q > 1");
   }

   SparseMatrix mat;

   if(q == 1)
   {
      SparseSM::Chebyshev::LinearMap::I1 qi(rows, cols, lb, ub);
      mat = qi.mat();
   }
   else if(q == 2)
   {
      SparseSM::Chebyshev::LinearMap::I2 qi(rows, cols, lb, ub);
      mat = qi.mat();
   }
   else if(q == 3)
   {
      SparseSM::Chebyshev::LinearMap::I3 qi(rows, cols, lb, ub);
      mat = qi.mat();
   }
   else if(q == 4)
   {
      SparseSM::Chebyshev::LinearMap::I4 qi(rows, cols, lb, ub);
      mat = qi.mat();
   }
   else
   {
      throw std::logic_error("Quasi-inverse of order > 4 is not implemented");
   }

   return mat;
}

SparseMatrix matIq(const int q, const int i, const int rows, const int cols, const Internal::MHDFloat lb, const Internal::MHDFloat ub)
{
   if(i > q)
   {
      throw std::logic_error("Cannot have negative quasi-inverse order: i > q");
   }

   SparseMatrix mat;
   if(q == i)
   {
         SparseSM::Chebyshev::LinearMap::Id qid(rows, cols, lb, ub, q);
         mat = qid.mat();
   }
   else if(i == 0)
   {
      mat = selectIq(q, rows, cols, lb, ub);
   }
   else
   {
      SparseSM::Chebyshev::LinearMap::Id qid(rows, rows, lb, ub, q);
      mat = qid.mat()*selectIq(q-i, rows, cols, lb, ub);
   }
   if(mat.rows() != rows || mat.cols() != cols)
   {
      throw std::logic_error("Someting whent wring");
   }

   return mat;
}

} // namespace Utils
} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
