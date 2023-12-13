/**
 * @file PyGeostrophic2Tor.cpp
 * @brief Source of the implementation of the projection operator from the geostrophic basis to Worland in Python
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <Eigen/Dense>

// Project includes
//
#include "PyGeostrophic2Tor.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/PyQuICC/CoreWrapper.hpp"
#include "QuICC/PyQuICC/Tools.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   PyGeostrophic2Tor::PyGeostrophic2Tor(const int nN, const int maxnl, const int nR, const int maxNug, const ArrayI& nli, const std::vector<int>& nIdx, const Scalar_t ugAlpha, const Scalar_t ugDBeta, const Scalar_t alpha, const Scalar_t dBeta, const int q)
      : IEmbeddedOperator(maxnl*nN, maxNug+1, alpha, dBeta), mNn(nN), mMaxnl(maxnl), mNr(nR), mMaxNug(maxNug), mNlist(nli), mNidx(nIdx)
   {
      PyQuICC::CoreWrapper::init();
      PyQuICC::CoreWrapper::import("quicc_solver.model.boussinesq.sphere.modifiedtaylor.linear.helper");
      PyQuICC::CoreWrapper::createClass("Helper");
   }

   PyGeostrophic2Tor::~PyGeostrophic2Tor()
   {
      PyQuICC::CoreWrapper::finalize();
   }

   void PyGeostrophic2Tor::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
   {
      switch(this->type())
      {
         case WorlandKind::CHEBYSHEV:
            this->buildChebyshevOp(mat, rows, cols);
            break;
         case WorlandKind::LEGENDRE:
            throw std::logic_error("Legendre basis operator not implemented");
            break;
         case WorlandKind::CYLENERGY:
            throw std::logic_error("Cylindrical energy basis operator not implemented");
            break;
         case WorlandKind::SPHENERGY:
            throw std::logic_error("Spherical energy basis operator not implemented");
            break;
      }
   }

   void PyGeostrophic2Tor::buildChebyshevOp(Internal::Matrix& mat, const int rows, const int cols) const
   {
      const bool isTriangular = true;
      const auto& nN = this->mNn;
      const auto& maxnl = this->mMaxnl;
      const auto& maxNug = this->mMaxNug;
      const auto& nIdx = this->mNidx;

      // Get simulation spectral resolution array
      ArrayI simRes(3);
      simRes << nN, maxnl, maxnl;

      Matrix pyMat = Matrix::Zero(nN*maxnl, maxNug+1);
      Matrix tmp(nN*maxnl, maxNug+1);
      for(int n: nIdx)
      {
         // Prepare Python call arguments
         PyObject *pArgs, *pValue;
         pArgs = PyTuple_New(2);

         // Get resolution
         pValue = PyQuICC::Tools::makeTuple(simRes);
         PyTuple_SetItem(pArgs, 0, pValue);
         PyTuple_SetItem(pArgs, 1, PyLong_FromLong(n));

         // Call model operator Python routine
         PyQuICC::CoreWrapper::setMethod("geo2tor_mat_n");
         pValue = PyQuICC::CoreWrapper::callMethod(pArgs);
         Py_DECREF(pArgs);

         // Convert Python matrix into Matrix
         PyQuICC::CoreWrapper::fillMatrix(tmp, pValue);
         pyMat.col(n) = tmp.col(n);
         Py_DECREF(pValue);

         // Finalise Python interpreter
         PyQuICC::CoreWrapper::cleanup();
      }

      mat = pyMat.cast<Internal::MHDFloat>();
   }

} // Worland
} // DenseSM
} // QuICC
