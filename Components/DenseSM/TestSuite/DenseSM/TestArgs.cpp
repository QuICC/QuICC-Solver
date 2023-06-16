/**
 * @file TestArgs.cpp
 * @brief Source of test arguments
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "TestSuite/DenseSM/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace DenseSM {

   TestArgs::TestArgs()
      : useDefault(true), dumpData(false), type(TestType::DENSE), ulp(11)
   {
   }

   void TestArgs::setType(const std::string& type)
   {
      if(type == "dense")
      {
         this->type = TestType::DENSE;
      }
      else if(type == "sparse")
      {
         this->type = TestType::SPARSE;
      }
      else
      {
         throw std::logic_error("Unsupported test type!");
      }
   }

}
}
}
