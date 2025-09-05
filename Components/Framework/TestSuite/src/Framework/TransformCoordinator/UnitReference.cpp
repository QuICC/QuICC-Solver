/**
 * @file UnitReference.cpp
 * @brief Generate unit spectrum reference input and ouput from file
 */

// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>
#include <fstream>
#include <limits>

// Project includes
//
#include "QuICC/TestSuite/Framework/TransformCoordinator/UnitReference.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace TCoord {

   UnitReference::UnitReference(const Test::FieldId id)
      : mFieldId(id)
   {
   }

   MHDComplex UnitReference::inScalar(Test& test, int i, int j, int k)
   {
      if(this->mFieldId == Test::FieldId::SCALAR ||
         this->mFieldId == Test::FieldId::SCALAR_AND_TORPOL)
      {
         return unitReference(test, i, j, k);
      }
      else
      {
         return MHDComplex(0.0);
      }
   }

   MHDComplex UnitReference::inTor(Test& test, int i, int j, int k)
   {
      if(this->mFieldId == Test::FieldId::TOR ||
         this->mFieldId == Test::FieldId::TORPOL ||
         this->mFieldId == Test::FieldId::SCALAR_AND_TORPOL)
      {
         return unitReference(test, i, j, k);
      }
      else
      {
         return MHDComplex(0.0);
      }
   }

   MHDComplex UnitReference::inPol(Test& test, int i, int j, int k)
   {
      if(this->mFieldId == Test::FieldId::POL ||
         this->mFieldId == Test::FieldId::TORPOL ||
         this->mFieldId == Test::FieldId::SCALAR_AND_TORPOL)
      {
         return unitReference(test, i, j, k);
      }
      else
      {
         return MHDComplex(0.0);
      }
   }

   MHDComplex UnitReference::refScalar(Test& test, int i, int j, int k)
   {
      if(this->mFieldId == Test::FieldId::SCALAR ||
         this->mFieldId == Test::FieldId::SCALAR_AND_TORPOL)
      {
         return unitReference(test, i, j, k);
      }
      else
      {
         return MHDComplex(0.0);
      }
   }

   MHDComplex UnitReference::refTor(Test& test, int i, int j, int k)
   {
      if(this->mFieldId == Test::FieldId::TOR ||
         this->mFieldId == Test::FieldId::TORPOL ||
         this->mFieldId == Test::FieldId::SCALAR_AND_TORPOL)
      {
         return unitReference(test, i, j, k);
      }
      else
      {
         return MHDComplex(0.0);
      }
   }

   MHDComplex UnitReference::refPol(Test& test, int i, int j, int k)
   {
      if(this->mFieldId == Test::FieldId::POL ||
         this->mFieldId == Test::FieldId::TORPOL ||
         this->mFieldId == Test::FieldId::SCALAR_AND_TORPOL)
      {
         return unitReference(test, i, j, k);
      }
      else
      {
         return MHDComplex(0.0);
      }
   }
}
}
}
}
