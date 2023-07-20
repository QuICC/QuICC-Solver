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
      : useDefault(true), dumpData(false), timeOnly(false), type(TestType::DENSE), ulp(11), np(0), rank(0), iter(1)
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
