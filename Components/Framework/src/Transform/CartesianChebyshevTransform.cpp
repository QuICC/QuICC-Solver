/** 
 * @file CartesianChebyshevTransform.cpp
 * @brief Source of the implementation of the transform
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/CartesianChebyshevTransform.hpp"
#include "QuICC/Transform/RegisterCartesianChebyshevMap.hpp"

namespace QuICC {

namespace Transform {

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
      for(const auto& f: RegisterCartesianChebyshevMap::mapper())
      {
         this->mImpl.addOperator(*f);
      }
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
