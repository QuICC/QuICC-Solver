/** 
 * @file AnnulusChebyshevTransform.cpp
 * @brief Source of the implementation of the annulus Chebyshev transform
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/AnnulusChebyshevTransform.hpp"

// Project includes
//

#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY2.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D.hpp"
//NOt IMPLEMENTED YET #include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY1D1.hpp"
//NOt IMPLEMENTED YET#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D1DivY1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"
                                              
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"

#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/Energy.hpp"

#include "QuICC/Transform/Forward/P.hpp"

#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/Overr1.hpp"
#include "QuICC/Transform/Backward/Overr2.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
//NOt IMPLEMENTED YET#include "QuICC/Transform/Backward/Overr1D1.hpp"
//NOt IMPLEMENTED YET#include "QuICC/Transform/Backward/D1Overr1.hpp"

#include "QuICC/Transform/Reductor/Energy.hpp"

namespace QuICC {

namespace Transform {

   AnnulusChebyshevTransform::AnnulusChebyshevTransform()
   {
   }

   AnnulusChebyshevTransform::~AnnulusChebyshevTransform()
   {
   }

   void AnnulusChebyshevTransform::requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const
   {
      this->mImpl.requiredOptions(list, dimId);
   }

   void AnnulusChebyshevTransform::setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId)
   {
      this->mImpl.setOptions(options, dimId);
   }

   Array AnnulusChebyshevTransform::meshGrid() const
   {
      return this->mImpl.meshGrid();
   }

   void AnnulusChebyshevTransform::init(AnnulusChebyshevTransform::SharedSetupType spSetup)
   {
      // Initialize transform implementation
      this->mImpl.init(spSetup);

      // Initialise operators
      this->initOperators();
   }

   void AnnulusChebyshevTransform::initOperators()
   {
      // Create projectors
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::P>(Backward::P::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::DivY1>(Backward::Overr1::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::DivY2>(Backward::Overr2::id());
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::D1>(Backward::D1::id());
      //NOT IMPLEMENTED YET this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::D1>(Backward::Overr1D1::id());
      //NOT IMPLEMENTED YET this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Projector::D1>(Backward::D1Overr1::id());

      // Create integrators
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Integrator::P>(Forward::P::id());

      // Create reductors
      this->mImpl.addOperator<Fft::Chebyshev::LinearMap::Reductor::Energy>(Reductor::Energy::id());
   }                                                 

   void AnnulusChebyshevTransform::forward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void AnnulusChebyshevTransform::backward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void AnnulusChebyshevTransform::reduce(Matrix& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   //
   // Disabled transforms
   //

   void AnnulusChebyshevTransform::forward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void AnnulusChebyshevTransform::forward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void AnnulusChebyshevTransform::forward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void AnnulusChebyshevTransform::backward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void AnnulusChebyshevTransform::backward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void AnnulusChebyshevTransform::backward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void AnnulusChebyshevTransform::reduce(MatrixZ&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void AnnulusChebyshevTransform::reduce(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void AnnulusChebyshevTransform::reduce(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   MHDFloat AnnulusChebyshevTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += this->mImpl.requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void AnnulusChebyshevTransform::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MHDFloat mem = this->mImpl.requiredStorage();

      StorageProfilerMacro_update(StorageProfilerMacro::TRAANNULUSCHEBYSHEV, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::TRANSFORMS, mem);
#endif // QUICC_STORAGEPROFILE
   }

}
}
