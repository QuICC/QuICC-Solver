/**
 * file SphereWorlandTransform.cpp
 * @brief Source of the implementation of the Worland transform in a sphere
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/SphereWorlandTransform.hpp"
#include "QuICC/Transform/RegisterSphereWorlandMap.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

   void SphereWorlandTransform::requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const
   {
      this->mImpl.requiredOptions(list, dimId);
   }

   void SphereWorlandTransform::setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId)
   {
      this->mImpl.setOptions(options, dimId);
   }

   Array SphereWorlandTransform::meshGrid() const
   {
      return this->mImpl.meshGrid();
   }

   void SphereWorlandTransform::init(SphereWorlandTransform::SharedSetupType spSetup)
   {
      // Initialize transform implementation
      this->mImpl.init(spSetup);

      // Initialise the ID to operator mapping
      this->initOperators();
   }

   void SphereWorlandTransform::initOperators()
   {
      for(const auto& f: RegisterSphereWorlandMap::mapper())
      {
         this->mImpl.addOperator(*f);
      }
   }

   void SphereWorlandTransform::forward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id)
   {
      Profiler::RegionFixture<3> fix("SphereWorlandTransform::forward");
      this->mImpl.transform(rOut, in, id);
   }

   void SphereWorlandTransform::reduce(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void SphereWorlandTransform::backward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id)
   {
      Profiler::RegionFixture<3> fix("SphereWorlandTransform::backward");
      this->mImpl.transform(rOut, in, id);
   }

   MHDFloat SphereWorlandTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += this->mImpl.requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void SphereWorlandTransform::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MHDFloat mem = this->mImpl.requiredStorage();

      StorageProfilerMacro_update(StorageProfilerMacro::TRASPHEREWORLAND, mem);
      StorageProfilerMacro_update(StorageProfilerMacro::TRANSFORMS, mem);
#endif // QUICC_STORAGEPROFILE
   }
} // Transform
} // QuICC
