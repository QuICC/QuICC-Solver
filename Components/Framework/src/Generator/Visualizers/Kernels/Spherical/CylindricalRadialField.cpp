/**
 * @file CylindricalRadialField.cpp
 * @brief Source of physical space kernel to compute cylindrical radial compoent of spherical field
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/Visualizers/Kernels/Spherical/CylindricalRadialField.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"
#include "QuICC/PhysicalOperators/SphericalSComponent.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Spherical {

   CylindricalRadialField::CylindricalRadialField()
      : IPhysicalKernel()
   {
   }

   CylindricalRadialField::~CylindricalRadialField()
   {
   }

   std::size_t CylindricalRadialField::name() const
   {
      return this->mName;
   }

   void CylindricalRadialField::setVector(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

      this->mName = name;

      this->setField(name, spField);
   }

   void CylindricalRadialField::init(FieldType::Id type, const MHDFloat scale)
   {
      this->mFieldType = type;

      this->mScale = scale;
   }

   void CylindricalRadialField::setMesh(std::shared_ptr<std::vector<Array> > spMesh)
   {
      IPhysicalKernel::setMesh(spMesh);

      this->mCosTheta = this->mspMesh->at(1).array().cos();
      this->mSinTheta = this->mspMesh->at(1).array().sin();
   }

   void CylindricalRadialField::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id) const
   {
      if(this->mFieldType == FieldType::VECTOR)
      {
         std::visit([&](auto&& p){Physical::SphericalSComponent::set(rNLComp, p->dom(0).res(), this->mCosTheta, this->mSinTheta, p->dom(0).phys(), this->mScale);}, this->vector(this->name()));
      } else if(this->mFieldType == FieldType::GRADIENT)
      {
         throw std::logic_error("Z Component of gradient not implemented yet!");

         //Physical::SphericalZComponent::add(rNLComp, this->res(), this->mCosTheta, this->mSinTheta, this->vector(mFieldName).dom(0).grad());
      } else if(this->mFieldType == FieldType::CURL)
      {
         std::visit([&](auto&& p){Physical::SphericalSComponent::set(rNLComp, p->dom(0).res(), this->mCosTheta, this->mSinTheta, p->dom(0).curl(), this->mScale);}, this->vector(this->name()));
      }
   }

}
}
}
}
