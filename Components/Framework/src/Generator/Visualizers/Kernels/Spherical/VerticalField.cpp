/**
 * @file VerticalField.cpp
 * @brief Source of physical space kernel to compute spherical vertical (Z) field component
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/Visualizers/Kernels/Spherical/VerticalField.hpp"

// Project includes
//
#include "Types/Math.hpp"
#include "QuICC/PhysicalNames/Undefined.hpp"
#include "QuICC/PhysicalOperators/SphericalZComponent.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Spherical {

   VerticalField::VerticalField()
      : IPhysicalKernel(), mNameId(PhysicalNames::Undefined::id()), mCompNameId(PhysicalNames::Undefined::id())
   {
   }

   VerticalField::~VerticalField()
   {
   }

   std::size_t VerticalField::nameId() const
   {
      return this->mNameId;
   }

   std::size_t VerticalField::compNameId() const
   {
      return this->mCompNameId;
   }

   bool VerticalField::hasScalar() const
   {
      return (this->mScalars.count(this->mCompNameId) == 1);
   }

   void VerticalField::setVector(std::size_t nameId, Framework::Selector::VariantSharedVectorVariable spField)
   {
      // Safety assertion
      assert(this->mVectors.count(nameId) == 0);

      this->mNameId = nameId;

      this->setField(nameId, spField);
   }

   void VerticalField::setScalar(std::size_t nameId, Framework::Selector::VariantSharedScalarVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(nameId) == 0);

      this->mCompNameId = nameId;

      this->setField(nameId, spField);
   }

   void VerticalField::init(FieldType::Id type, const MHDFloat scale)
   {
      this->mFieldType = type;

      this->mScale = scale;
   }

   void VerticalField::setMesh(std::shared_ptr<std::vector<Array> > spMesh)
   {
      IPhysicalKernel::setMesh(spMesh);

      this->mCosTheta = this->mspMesh->at(1).array().cos();
      this->mSinTheta = this->mspMesh->at(1).array().sin();
   }

   void VerticalField::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id) const
   {
      if(this->mFieldType == FieldType::VECTOR)
      {
         if(this->hasScalar())
         {
            std::visit([&](auto&& p, auto&& q){Physical::SphericalZComponent::set(q->rDom(0).rPhys(), p->dom(0).res(), this->mCosTheta, this->mSinTheta, p->dom(0).phys(), this->mScale);}, this->vector(this->nameId()), this->scalar(this->compNameId()));
         }
         else
         {
            std::visit([&](auto&& p){Physical::SphericalZComponent::set(rNLComp, p->dom(0).res(), this->mCosTheta, this->mSinTheta, p->dom(0).phys(), this->mScale);}, this->vector(this->nameId()));
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
            std::visit([&](auto&& p, auto&& q){Physical::SphericalZComponent::set(q->rDom(0).rPhys(), p->dom(0).res(), this->mCosTheta, this->mSinTheta, p->dom(0).curl(), this->mScale);}, this->vector(this->nameId()), this->scalar(this->compNameId()));
         }
         else
         {
            std::visit([&](auto&& p){Physical::SphericalZComponent::set(rNLComp, p->dom(0).res(), this->mCosTheta, this->mSinTheta, p->dom(0).curl(), this->mScale);}, this->vector(this->nameId()));
         }
      }
   }

}
}
}
}
