/**
 * @file ChebyshevEnergy.cpp
 * @brief Source of the interface for a generic API for Chebyshev energy reductor
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Backend/ChebyshevEnergy.hpp"
#if defined QUICC_FFT_CHEBYSHEV_FFTW
   #include "QuICC/Transform/Fft/Backend/Fftw/ChebyshevEnergy.hpp"
   #define BACKENDIMPL Fftw
#elif defined QUICC_FFT_CHEBYSHEV_CUFFT
   #include "QuICC/Transform/Fft/Backend/CuFft/ChebyshevEnergy.hpp"
   #define BACKENDIMPL CuFft
#endif

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   struct ChebyshevEnergy::BackendImpl: public BACKENDIMPL::ChebyshevEnergy
   {
      BackendImpl() = default;
   };

   ChebyshevEnergy::ChebyshevEnergy()
   {
      this->mpImpl = std::make_shared<BackendImpl>();
   }

   void ChebyshevEnergy::init(const SetupType& setup) const
   {
      this->mpImpl->init(setup);
   }

   void ChebyshevEnergy::setScaler(const Array& scaler) const
   {
      this->mpImpl->setScaler(scaler);
   }
   
   // anelastic overload
   void ChebyshevEnergy::init(const SetupType& setup, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const
   {
      this->mpImpl->init(setup, pF);
   }

   void ChebyshevEnergy::setSpectralOperator(const SparseMatrix& mat) const
   {
      this->mpImpl->setSpectralOperator(mat);
   }

   void ChebyshevEnergy::input(Matrix& tmp, const Matrix& in) const
   {
      this->mpImpl->input(tmp, in);
   }

   void ChebyshevEnergy::input(Matrix& tmp, const Matrix& in,
      const int shift) const
   {
      this->mpImpl->input(tmp, in, shift);
   }

   void ChebyshevEnergy::input(Matrix& tmp, const MatrixZ& in,
      const bool useReal) const
   {
      this->mpImpl->input(tmp, in, useReal);
   }

   void ChebyshevEnergy::input(Matrix& tmp, const MatrixZ& in,
      const int shift, const bool useReal) const
   {
      this->mpImpl->input(tmp, in, shift, useReal);
   }

   void ChebyshevEnergy::square(Matrix& tmp, const Matrix& in,const bool isFirst) const
   {
      this->mpImpl->square(tmp, in, isFirst);
   }

   void ChebyshevEnergy::output(Eigen::Ref<Matrix> rOut, const Matrix& tmp) const
   {
      this->mpImpl->output(rOut, tmp);
   }

   void ChebyshevEnergy::outputGrid(Eigen::Ref<Matrix> rOut, const Matrix& tmp) const
   {
      this->mpImpl->outputGrid(rOut, tmp);
   }

   void ChebyshevEnergy::outputSpectral(Eigen::Ref<Matrix> rOut, const Matrix& tmp) const
   {
      this->mpImpl->outputSpectral(rOut, tmp);
   }

   void ChebyshevEnergy::applyFft(Eigen::Ref<Matrix> phys, const Eigen::Ref<const Matrix>& mods) const
   {
      this->mpImpl->applyFft(phys, mods);
   }

   void ChebyshevEnergy::applyFwdFft(Eigen::Ref<Matrix> mods, const Eigen::Ref<const Matrix>& phys) const
   {
      this->mpImpl->applyFwdFft(mods, phys);
   }

   void ChebyshevEnergy::addSolver(const int extraRows) const
   {
      this->mpImpl->addSolver(extraRows);
   }

   void ChebyshevEnergy::getSolution(Matrix& tmp, const int zeroRows, const int extraRows) const
   {
      this->mpImpl->getSolution(tmp, zeroRows, extraRows);
   }

   Matrix& ChebyshevEnergy::getStorage(const StorageKind kind) const
   {
      return this->mpImpl->getStorage(kind);
   }

   Fftw::DifferentialSolver& ChebyshevEnergy::solver() const
   {
      return this->mpImpl->solver();
   }

   Array& ChebyshevEnergy::getEGrid() const
   {
      return this->mpImpl->getEGrid();
   }

   MHDFloat ChebyshevEnergy::getFftScaling() const
   {
      return this->mpImpl->getFftScaling();
   }

   void ChebyshevEnergy::setExtraSize(int extraSize) const
   {
      this->mpImpl->setExtraSize(extraSize);
   }

   int ChebyshevEnergy::getExtraSize() const
   {
      return this->mpImpl->getExtraSize();
   }
}
}
}
}
