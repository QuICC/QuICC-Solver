/**
 * file SphereBesselTransform.cpp
 * @brief Source of the implementation of the spherical Bessel transform in a sphere
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/SphereBesselTransform.hpp"
#include "QuICC/Transform/RegisterSphereBesselMap.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

   void SphereBesselTransform::requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const
   {
      this->mImpl.requiredOptions(list, dimId);
   }

   void SphereBesselTransform::setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId)
   {
      this->mImpl.setOptions(options, dimId);
   }

   Array SphereBesselTransform::meshGrid() const
   {
      return this->mImpl.meshGrid();
   }

   void SphereBesselTransform::init(SphereBesselTransform::SharedSetupType spSetup)
   {
      // Initialize transform implementation
      this->mImpl.init(spSetup);

      // Initialise the ID to operator mapping
      this->initOperators();
   }

   void SphereBesselTransform::initOperators()
   {
      for(const auto& f: RegisterSphereBesselMap::mapper())
      {
         this->mImpl.addOperator(*f);
      }
   }

   void SphereBesselTransform::forward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      Profiler::RegionFixture<3> fix("SphereBesselTransform::forward");
      this->mImpl.transform(rOut, in, id);
   }

   void SphereBesselTransform::reduce(Matrix& rOut, const MatrixZ& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void SphereBesselTransform::backward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id)
   {
      Profiler::RegionFixture<3> fix("SphereBesselTransform::backward");
      this->mImpl.transform(rOut, in, id);
   }

   //
   // Disabled transforms
   //

   void SphereBesselTransform::forward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereBesselTransform::forward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereBesselTransform::forward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereBesselTransform::backward(Matrix&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereBesselTransform::backward(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereBesselTransform::backward(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereBesselTransform::reduce(MatrixZ&, const MatrixZ&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereBesselTransform::reduce(MatrixZ&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   void SphereBesselTransform::reduce(Matrix&, const Matrix&, const std::size_t)
   {
      this->unimplemented();
   }

   MHDFloat SphereBesselTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += this->mImpl.requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void SphereBesselTransform::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MHDFloat mem = this->mImpl.requiredStorage();

      StorageProfilerMacro_update(StorageProfilerMacro::TRASPHEREBESSEL, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::TRANSFORMS, mem);
#endif // QUICC_STORAGEPROFILE
   }
} // Transform
} // QuICC
