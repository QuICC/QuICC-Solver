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
#include "QuICC/PhysicalNames/Undefined.hpp"
#include "QuICC/PhysicalOperators/SphericalSComponent.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Spherical {

   CylindricalRadialField::CylindricalRadialField()
      : IPhysicalKernel(), mNameId(PhysicalNames::Undefined::id()), mCompNameId(PhysicalNames::Undefined::id())
   {
   }

   CylindricalRadialField::~CylindricalRadialField()
   {
   }

   std::size_t CylindricalRadialField::nameId() const
   {
      return this->mNameId;
   }

   std::size_t CylindricalRadialField::compNameId() const
   {
      return this->mCompNameId;
   }

   bool CylindricalRadialField::hasScalar() const
   {
      return (this->mScalars.count(this->mCompNameId) == 1);
   }

   void CylindricalRadialField::setVector(std::size_t nameId, Framework::Selector::VariantSharedVectorVariable spField)
   {
      // Safety assertion
      assert(this->mVectors.count(nameId) == 0);

      this->mNameId = nameId;

      this->setField(nameId, spField);
   }

   void CylindricalRadialField::setScalar(std::size_t nameId, Framework::Selector::VariantSharedScalarVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(nameId) == 0);

      this->mCompNameId = nameId;

      this->setField(nameId, spField);
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
         if(this->hasScalar())
         {
            std::visit([&](auto&& p, auto&& q){Physical::SphericalSComponent::set(q->rDom(0).rPhys(), p->dom(0).res(), this->mCosTheta, this->mSinTheta, p->dom(0).phys(), this->mScale);}, this->vector(this->nameId()), this->scalar(this->compNameId()));
         }
         else
         {
            std::visit([&](auto&& p){Physical::SphericalSComponent::set(rNLComp, p->dom(0).res(), this->mCosTheta, this->mSinTheta, p->dom(0).phys(), this->mScale);}, this->vector(this->nameId()));
         }
      }
      else if(this->mFieldType == FieldType::GRADIENT)
      {
         throw std::logic_error("Z Component of gradient not implemented yet!");

      }
      else if(this->mFieldType == FieldType::CURL)
      {
         if(this->hasScalar())
         {
            std::visit([&](auto&& p, auto&& q){Physical::SphericalSComponent::set(q->rDom(0).rPhys(), p->dom(0).res(), this->mCosTheta, this->mSinTheta, p->dom(0).curl(), this->mScale);}, this->vector(this->nameId()), this->scalar(this->compNameId()));
         }
         else
         {
            std::visit([&](auto&& p){Physical::SphericalSComponent::set(rNLComp, p->dom(0).res(), this->mCosTheta, this->mSinTheta, p->dom(0).curl(), this->mScale);}, this->vector(this->nameId()));
         }
      }
   }

}
}
}
}
