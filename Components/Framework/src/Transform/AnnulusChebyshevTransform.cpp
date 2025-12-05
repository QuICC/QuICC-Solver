/** 
 * @file AnnulusChebyshevTransform.cpp
 * @brief Source of the implementation of the annulus Chebyshev transform
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/AnnulusChebyshevTransform.hpp"
#include "QuICC/Transform/RegisterAnnulusChebyshevMap.hpp"

namespace QuICC {

namespace Transform {

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
      for(const auto& f: RegisterAnnulusChebyshevMap::mapper())
      {
         this->mImpl.addOperator(*f);
      }
   }                                                 

   void AnnulusChebyshevTransform::forward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void AnnulusChebyshevTransform::backward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
   }

   void AnnulusChebyshevTransform::reduce(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id)
   {
      this->mImpl.transform(rOut, in, id);
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
