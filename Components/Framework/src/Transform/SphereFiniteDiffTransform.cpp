/**
 * file SphereFiniteDiffTransform.cpp
 * @brief Source of the implementation of the Finite Differences transform in a sphere
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/SphereFiniteDiffTransform.hpp"
#include "QuICC/Transform/RegisterSphereFiniteDiffMap.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

   void SphereFiniteDiffTransform::requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const
   {
      this->mImpl.requiredOptions(list, dimId);
   }

   void SphereFiniteDiffTransform::setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId)
   {
      this->mImpl.setOptions(options, dimId);
   }

   Array SphereFiniteDiffTransform::meshGrid() const
   {
      return this->mImpl.meshGrid();
   }

   void SphereFiniteDiffTransform::init(SphereFiniteDiffTransform::SharedSetupType spSetup)
   {
      // Initialize transform implementation
      this->mImpl.init(spSetup);

      // Initialise the ID to operator mapping
      this->initOperators();
   }

   void SphereFiniteDiffTransform::initOperators()
   {
      for(const auto& f: RegisterSphereFiniteDiffMap::mapper())
      {
         this->mImpl.addOperator(*f);
      }
   }

   void SphereFiniteDiffTransform::forward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      Profiler::RegionFixture<3> fix("SphereFiniteDiffTransform::forward");
      this->mImpl.transform(rOut, in, id);
   }

   void SphereFiniteDiffTransform::reduce(Matrix& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void SphereFiniteDiffTransform::backward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      Profiler::RegionFixture<3> fix("SphereFiniteDiffTransform::backward");
      this->mImpl.transform(rOut, in, id);
   }

   //
   // Disabled transforms
   //

   void SphereFiniteDiffTransform::forward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereFiniteDiffTransform::forward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereFiniteDiffTransform::forward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereFiniteDiffTransform::backward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereFiniteDiffTransform::backward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereFiniteDiffTransform::backward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereFiniteDiffTransform::reduce(MatrixZ&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereFiniteDiffTransform::reduce(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereFiniteDiffTransform::reduce(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   MHDFloat SphereFiniteDiffTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += this->mImpl.requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void SphereFiniteDiffTransform::profileStorage() const
   {
   }
} // Transform
} // QuICC
