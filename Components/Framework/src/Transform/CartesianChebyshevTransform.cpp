/** 
 * @file CartesianChebyshevTransform.cpp
 * @brief Source of the implementation of the transform
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/CartesianChebyshevTransform.hpp"

// Project includes
//

#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D4.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D3.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D2.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"

#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I4.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2D1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I4D1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2_I2D1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I4D1_I2.hpp"

#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/Energy.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/EnergyD1.hpp"

#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D2.hpp"

#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/I2P.hpp"
#include "QuICC/Transform/Forward/I2D1.hpp"
#include "QuICC/Transform/Forward/I2ZI2D1.hpp"
#include "QuICC/Transform/Forward/I4P.hpp"
#include "QuICC/Transform/Forward/I4D1.hpp"
#include "QuICC/Transform/Forward/I4D1ZI2.hpp"

#include "QuICC/Transform/Reductor/Energy.hpp"
#include "QuICC/Transform/Reductor/EnergyD1.hpp"

namespace QuICC {

namespace Transform {

   CartesianChebyshevTransform::CartesianChebyshevTransform()
   {
   }

   CartesianChebyshevTransform::~CartesianChebyshevTransform()
   {
   }

   void CartesianChebyshevTransform::requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const
   {
      this->mImpl.requiredOptions(list, dimId);
   }

   void CartesianChebyshevTransform::setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId)
   {
      this->mImpl.setOptions(options, dimId);
   }

   Array CartesianChebyshevTransform::meshGrid() const
   {
      return this->mImpl.meshGrid();
   }

   void CartesianChebyshevTransform::init(CartesianChebyshevTransform::SharedSetupType spSetup)
   {
      // Initialize transform implementation
      this->mImpl.init(spSetup);

      // Initialise FFTW plans
      this->initOperators();
   }

   void CartesianChebyshevTransform::initOperators()
   {
      // Create projectors
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::D1>(Backward::D1::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::D2>(Backward::D2::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::P>(Backward::P::id());

      // Create integrators
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Integrator::P>(Forward::P::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Integrator::I2>(Forward::I2P::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Integrator::I4>(Forward::I4P::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Integrator::I2D1>(Forward::I2D1::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Integrator::I4D1>(Forward::I4D1::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Integrator::I2_I2D1>(Forward::I2ZI2D1::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Integrator::I4D1_I2>(Forward::I4D1ZI2::id());

      // Create reductors
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Reductor::Energy>(Reductor::Energy::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Reductor::EnergyD1>(Reductor::EnergyD1::id());
   }                                                 

   void CartesianChebyshevTransform::forward(Matrix& rOut, const Matrix& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void CartesianChebyshevTransform::backward(Matrix& rOut, const Matrix& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void CartesianChebyshevTransform::forward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void CartesianChebyshevTransform::backward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void CartesianChebyshevTransform::reduce(Matrix& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void CartesianChebyshevTransform::reduce(Matrix& rOut, const Matrix& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   //
   // Disabled transforms
   //

   void CartesianChebyshevTransform::forward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void CartesianChebyshevTransform::forward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void CartesianChebyshevTransform::backward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void CartesianChebyshevTransform::backward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void CartesianChebyshevTransform::reduce(MatrixZ&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void CartesianChebyshevTransform::reduce(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   MHDFloat CartesianChebyshevTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += this->mImpl.requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void CartesianChebyshevTransform::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MHDFloat mem = this->mImpl.requiredStorage();

      StorageProfilerMacro_update(StorageProfilerMacro::TRACARTESIANCHEBYSHEV, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::TRANSFORMS, mem);
#endif // QUICC_STORAGEPROFILE
   }

}
}
