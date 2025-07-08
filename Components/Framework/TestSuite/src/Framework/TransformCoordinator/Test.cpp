/**
 * @file Test.cpp
 * @brief High level test setup
 */

// System includes
//
#include <catch2/catch.hpp>
#include <fstream>
#include <limits>

// Project includes
//
#include "QuICC/TestSuite/Framework/TransformCoordinator/Test.hpp"
#include "QuICC/TestSuite/Framework/TransformCoordinator/UnitReference.hpp"
#include "QuICC/TestSuite/Framework/TransformCoordinator/FileReference.hpp"
#include "TestSuite/Io.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace TCoord {

   Test::Test()
      : epsilon(std::numeric_limits<MHDFloat>::epsilon()), maxUlp(11)
   {
   }

   MHDFloat Test::tolerance() const
   {
      return this->maxUlp*this->epsilon;
   }

   void Test::configure(const int id, const std::string subdir)
   {
      const std::string refdir =  "_refdata/Framework/TransformCoordinator/" + subdir + "/";

      // Static setups
      if(id < 10)
      {
         switch(id)
         {
            case 0:
               this->fieldId = FieldId::SCALAR;
               this->kernelId = KernelId::PASSTHROUGH;
               this->pathId = PathId::BFLOOP;
               this->spRef = std::make_shared<UnitReference>(this->fieldId);
               break;
            case 1:
               this->fieldId = FieldId::TOR;
               this->kernelId = KernelId::PASSTHROUGH;
               this->pathId = PathId::BFLOOP;
               this->spRef = std::make_shared<UnitReference>(this->fieldId);
               break;
            case 2:
               this->fieldId = FieldId::POL;
               this->kernelId = KernelId::PASSTHROUGH;
               this->pathId = PathId::BFLOOP;
               this->spRef = std::make_shared<UnitReference>(this->fieldId);
               break;
            case 3:
               this->fieldId = FieldId::TORPOL;
               this->kernelId = KernelId::PASSTHROUGH;
               this->pathId = PathId::BFLOOP;
               this->spRef = std::make_shared<UnitReference>(this->fieldId);
               break;
            case 4:
               this->fieldId = FieldId::SCALAR_AND_TORPOL;
               this->kernelId = KernelId::PASSTHROUGH;
               this->pathId = PathId::BFLOOP;
               this->spRef = std::make_shared<UnitReference>(this->fieldId);
               break;
            default:
               throw std::logic_error("Undefined test case was requested");
         }
      }
      // Metadata file based setup
      else
      {
         this->fbase = refdir + "transform_loop_id" + std::to_string(id);
         this->spRef = std::make_shared<FileReference>(this->fbase);

         Array meta;
         std::string path = this->fbase + "_meta.dat";
         readList(meta, path);

         // Set field ID
         this->fieldId = static_cast<Test::FieldId>(static_cast<int>(meta(3)));
         // Set kernel ID
         this->kernelId = static_cast<Test::KernelId>(static_cast<int>(meta(4)));
         // Set path ID
         this->pathId = static_cast<Test::PathId>(static_cast<int>(meta(5)));
      }
   }
}
}
}
}
