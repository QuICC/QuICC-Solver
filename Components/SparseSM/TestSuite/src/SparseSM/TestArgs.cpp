/**
 * @file TestArgs.cpp
 * @brief Source of test arguments
 */

// Configuration includes
//

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
      else
      {
         throw std::logic_error("Unsupported test type!");
      }
   }

}
}
}
