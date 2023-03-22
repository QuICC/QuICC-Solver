/**
 * @file TestArgs.cpp
 * @brief Source of test arguments
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/TestSuite/SparseSM/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace SparseSM {

   TestArgs::TestArgs()
      : useDefault(true), dumpData(false), type(TestType::SPARSE), ulp(11)
   {
   }

   void TestArgs::setType(const std::string& type)
   {
      if(type == "sparse")
      {
         this->type = TestType::SPARSE;
      }
      else if(type == "banded")
      {
         this->type = TestType::BANDED;
      }
      else if(type == "boundary")
      {
         this->type = TestType::BOUNDARY;
      }
      else if(type == "stencil")
      {
         this->type = TestType::STENCIL;
      }
      else
      {
         throw std::logic_error("Unsupported test type!");
      }
   }

}
}
}
