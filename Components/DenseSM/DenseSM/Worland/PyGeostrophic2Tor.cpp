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
#include "DenseSM/Worland/details/GeostrophicTools.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/PyQuICC/CoreWrapper.hpp"
#include "QuICC/PyQuICC/Tools.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   PyGeostrophic2Tor::PyGeostrophic2Tor(const int nN, const int nL, const std::vector<int>& nIdx, const Scalar_t ugAlpha, const Scalar_t ugBeta, const bool isGenericBasis, const Scalar_t alpha, const Scalar_t dBeta, const bool isTriangular)
      : IWorlandOperator(nL*nN, nN, alpha, dBeta), mNn(nN), mNl(nL), mNidx(nIdx), mIsTriangular(isTriangular)
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
      const auto& nN = this->mNn;
      const auto& nL = this->mNl;
      const auto& nR = DenseSM::Worland::details::GeostrophicTools::cylTruncNr(nL, this->mIsTriangular);
      int nNug;
      if(this->mIsTriangular)
      {
         nNug = details::GeostrophicTools::cylTruncNug(nL, this->mIsTriangular);
      }
      else
      {
         nNug = details::GeostrophicTools::cylTruncNugC(nN, nL);
      }
      const auto nli = details::GeostrophicTools::nlist(nNug - 1, nL);
      const auto& nIdx = this->mNidx;

      // Get simulation spectral resolution array
      ArrayI simRes(3);
      simRes << nN, nL, nL;

      Matrix pyMat = Matrix::Zero(nN*nL, nNug);
      Matrix tmp(nN*nL, nNug);
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
