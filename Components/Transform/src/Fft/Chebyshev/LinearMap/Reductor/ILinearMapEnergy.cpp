/**
 * @file ILinearMapEnergy.cpp
 * @brief Source of the interface for a generic FFT based Chebyshev energy
 * reductor
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/ILinearMapEnergy.hpp"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

void ILinearMapEnergy::initBackend() const
{
   this->mBackend.init(*this->mspSetup);
}

// anelastic version (overload not possible without having to alter a bunch of other classes)
void ILinearMapEnergy::initBackendAnelastic(std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const
{
   // Set the extra size BEFORE backend initialization
   this->mExtraSize=pF->nN();
   // set it in the backend too
   this->mBackend.setExtraSize(pF->nN()); 

   this->mBackend.init(*this->mspSetup, pF);

}

void ILinearMapEnergy::transform(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in) const
{
   assert(this->isInitialized());
   assert(rOut.cols() == this->outCols());
   assert(rOut.rows() == this->outRows());

   auto& tmpIn = this->mBackend.getStorage(StorageKind::in);
   auto& tmpOut = this->mBackend.getStorage(StorageKind::out);
   auto& tmpSquare = this->mBackend.getStorage(StorageKind::mid);
   this->applyPreOperator(tmpIn, in, true);
   this->mBackend.applyFft(tmpOut, tmpIn);
   this->mBackend.square(tmpSquare, tmpOut, true);
   this->applyPreOperator(tmpIn, in, false);
   this->mBackend.applyFft(tmpOut, tmpIn);
   this->mBackend.square(tmpSquare, tmpOut, false);
   this->mBackend.applyFwdFft(tmpOut, tmpSquare);
   this->applyPostOperator(rOut, tmpOut);
}

// anelastic version:
void ILinearMapEnergy::transform(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const
{
   assert(this->isInitialized());
   assert(rOut.cols() == this->outCols());
   assert(rOut.rows() == this->outRows());

   auto eGrid = this->mBackend.getEGrid();
   
   auto rho = pF->evaluateLP(eGrid,0,0);

   auto& tmpIn = this->mBackend.getStorage(StorageKind::in); 
   auto& tmpOut = this->mBackend.getStorage(StorageKind::out); 
   auto& tmpSquare = this->mBackend.getStorage(StorageKind::mid);
   

   this->applyPreOperator(tmpIn, in, true); 
   this->mBackend.applyFft(tmpOut, tmpIn);
   this->mBackend.square(tmpSquare, tmpOut, true);

   this->applyPreOperator(tmpIn, in, false);
   this->mBackend.applyFft(tmpOut, tmpIn); 
   this->mBackend.square(tmpSquare, tmpOut, false); 
   
   tmpSquare = tmpSquare.array().colwise() / rho.array(); // divides energy by rho
   this->mBackend.applyFwdFft(tmpOut, tmpSquare);
   this->applyPostOperator(rOut, tmpOut);
}



void ILinearMapEnergy::transform(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const Matrix>& in) const
{
   assert(this->isInitialized());
   assert(rOut.cols() == this->outCols());
   assert(rOut.rows() == this->outRows());

   auto& tmpIn = this->mBackend.getStorage(StorageKind::in);
   auto& tmpOut = this->mBackend.getStorage(StorageKind::out);
   auto& tmpSquare = this->mBackend.getStorage(StorageKind::mid);
   this->applyPreOperator(tmpIn, in);
   this->mBackend.applyFft(tmpOut, tmpIn);
   this->mBackend.square(tmpSquare, tmpOut, true);
   this->mBackend.applyFwdFft(tmpOut, tmpSquare);
   this->applyPostOperator(rOut, tmpOut);
}

void ILinearMapEnergy::transform(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const Matrix>& in, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const
{
   assert(this->isInitialized());
   assert(rOut.cols() == this->outCols());
   assert(rOut.rows() == this->outRows());

   auto eGrid = this->mBackend.getEGrid();
   
   auto rho = pF->evaluateLP(eGrid,0,0);

   auto& tmpIn = this->mBackend.getStorage(StorageKind::in);
   auto& tmpOut = this->mBackend.getStorage(StorageKind::out);
   auto& tmpSquare = this->mBackend.getStorage(StorageKind::mid);
   this->applyPreOperator(tmpIn, in);
   this->mBackend.applyFft(tmpOut, tmpIn);
   this->mBackend.square(tmpSquare, tmpOut, true);
   tmpSquare = tmpSquare.array().colwise() / rho.array(); // divides energy by rho
   this->mBackend.applyFwdFft(tmpOut, tmpSquare);
   this->applyPostOperator(rOut, tmpOut);
}

void ILinearMapEnergy::transform(Eigen::Ref<MatrixZ>, const Eigen::Ref<const MatrixZ>&) const
{
   throw std::logic_error(
      "Data is not compatible with Chebyshev FFT energy reductor");
}

void ILinearMapEnergy::transform(Eigen::Ref<MatrixZ>, const Eigen::Ref<const Matrix>&) const
{
   throw std::logic_error(
      "Data is not compatible with Chebyshev FFT energy reductor");
}

MHDFloat ILinearMapEnergy::requiredStorage() const
{
   MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
   mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

   return mem;
}

int ILinearMapEnergy::outRows() const
{
   return this->mspSetup->blockSize();
}

int ILinearMapEnergy::outCols() const
{
   return 1;
}

} // namespace Reductor
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC
