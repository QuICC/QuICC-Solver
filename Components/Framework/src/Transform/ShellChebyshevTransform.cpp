/** 
 * @file ShellChebyshevTransform.cpp
 * @brief Source of the implementation of the spherical shell Chebyshev transform
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/ShellChebyshevTransform.hpp"

// Project includes
//

#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY2.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D1Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY1D1Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/SphRadLapl.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"
                                              
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I4Y3_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I4Y3D1Y1_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2Y2_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2Y1_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2Y1D1Y1_Zero.hpp"

#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/Energy.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/EnergyY2.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/EnergyD1Y1.hpp"

#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/R1.hpp"
#include "QuICC/Transform/Forward/Pol.hpp"
#include "QuICC/Transform/Forward/I4Q.hpp"
#include "QuICC/Transform/Forward/I4S.hpp"
#include "QuICC/Transform/Forward/I2Q.hpp"
#include "QuICC/Transform/Forward/I2S.hpp"
#include "QuICC/Transform/Forward/I2T.hpp"

#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/Overr1.hpp"
#include "QuICC/Transform/Backward/Overr2.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D1R1.hpp"
#include "QuICC/Transform/Backward/D2.hpp"
#include "QuICC/Transform/Backward/Overr1D1R1.hpp"
#include "QuICC/Transform/Backward/Slaplr.hpp"

#include "QuICC/Transform/Reductor/Energy.hpp"
#include "QuICC/Transform/Reductor/EnergyR2.hpp"
#include "QuICC/Transform/Reductor/EnergyD1R1.hpp"

namespace QuICC {

namespace Transform {

   ShellChebyshevTransform::ShellChebyshevTransform()
   {
   }

   ShellChebyshevTransform::~ShellChebyshevTransform()
   {
   }

   void ShellChebyshevTransform::requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const
   {
      this->mImpl.requiredOptions(list, dimId);
   }

   void ShellChebyshevTransform::setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId)
   {
      this->mImpl.setOptions(options, dimId);
   }

   Array ShellChebyshevTransform::meshGrid() const
   {
      return this->mImpl.meshGrid();
   }

   void ShellChebyshevTransform::init(ShellChebyshevTransform::SharedSetupType spSetup)
   {
      // Initialize transform implementation
      this->mImpl.init(spSetup);

      // Initialise operator
      this->initOperators();
   }

   void ShellChebyshevTransform::initOperators()
   {
      // Create projectors
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::P>(Backward::P::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::DivY1>(Backward::Overr1::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::DivY2>(Backward::Overr2::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::D1>(Backward::D1::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::D1Y1>(Backward::D1R1::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::D2>(Backward::D2::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::DivY1D1Y1>(Backward::Overr1D1R1::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::SphRadLapl>(Backward::Slaplr::id());

      // Create integrators
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Integrator::P>(Forward::P::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Integrator::Y1>(Forward::Pol::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Integrator::I4Y3_Zero>(Forward::I4Q::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Integrator::I4Y3D1Y1_Zero>(Forward::I4S::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Integrator::I2Y2_Zero>(Forward::I2T::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Integrator::I2Y1_Zero>(Forward::I2Q::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Integrator::I2Y1D1Y1_Zero>(Forward::I2S::id());

      // Create reductors
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Reductor::Energy>(Reductor::Energy::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Reductor::EnergyD1Y1>(Reductor::EnergyD1R1::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Reductor::EnergyY2>(Reductor::EnergyR2::id());
   }                                                 

   void ShellChebyshevTransform::forward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void ShellChebyshevTransform::backward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void ShellChebyshevTransform::reduce(Matrix& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   //
   // Disabled transforms
   //

   void ShellChebyshevTransform::forward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void ShellChebyshevTransform::forward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void ShellChebyshevTransform::forward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void ShellChebyshevTransform::backward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void ShellChebyshevTransform::backward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void ShellChebyshevTransform::backward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void ShellChebyshevTransform::reduce(MatrixZ&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void ShellChebyshevTransform::reduce(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void ShellChebyshevTransform::reduce(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   MHDFloat ShellChebyshevTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += this->mImpl.requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void ShellChebyshevTransform::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MHDFloat mem = this->mImpl.requiredStorage();

      StorageProfilerMacro_update(StorageProfilerMacro::TRASHELLCHEBYSHEV, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::TRANSFORMS, mem);
#endif // QUICC_STORAGEPROFILE
   }
}
}
