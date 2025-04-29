/**
 * @file IPfSolveBackend.cpp
 * @brief Source of the interface for a generic PFSOLVE based backend
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "QuICC/Transform/Fft/Backend/PfSolve/IPfSolveBackend.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace PfSolve {

IPfSolveBackend::IPfSolveBackend() : mPlan(NULL) { this->initLibrary(); }

IPfSolveBackend::~IPfSolveBackend() { this->cleanupFft(); }

void IPfSolveBackend::initLibrary() const {
  // PFSOLVE Fixture
  QuICC::Fft::Fftw::Library::getInstance();
}

void IPfSolveBackend::applyFft(Matrix &phys, const Matrix &mods) const {
  fftw_execute_r2r(this->mPlan, const_cast<MHDFloat *>(mods.data()),
                   phys.data());
}

void IPfSolveBackend::cleanupFft() {
  // Destroy plan
  if (this->mPlan) {
    fftw_destroy_plan(this->mPlan);
  }
}

MHDFloat IPfSolveBackend::requiredStorage() const {
  MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
  mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

  return mem;
}

void IPfSolveBackend::applyFft() const {
  std::logic_error("Backend not implemented.");
};
void IPfSolveBackend::applyFft(Matrix &, const MatrixZ &) const {
  std::logic_error("Backend not implemented.");
};
void IPfSolveBackend::applyFft(MatrixZ &, const Matrix &) const {
  std::logic_error("Backend not implemented.");
};
void IPfSolveBackend::applyFft(MatrixZ &, const MatrixZ &) const {
  std::logic_error("Backend not implemented.");
};

} // namespace PfSolve
} // namespace Backend
} // namespace Fft
} // namespace Transform
} // namespace QuICC
