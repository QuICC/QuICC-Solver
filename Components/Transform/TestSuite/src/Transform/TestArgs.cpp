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
#include "QuICC/TestSuite/Transform/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Transform {

   TestArgs::TestArgs()
      : useDefault(true), dumpData(false), type(TestType::PROJECTOR), ulp(11), np(0), rank(0)
   {
   }

   void TestArgs::setType(const std::string& type)
   {
      if(type == "projector")
      {
         this->type = TestType::PROJECTOR;
      }
      else if(type == "integrator")
      {
         this->type = TestType::INTEGRATOR;
      }
      else if(type == "reductor")
      {
         this->type = TestType::REDUCTOR;
      }
      else if(type == "bfloop")
      {
         this->type = TestType::BFLOOP;
      }
      else
      {
         throw std::logic_error("Unsupported test type!");
      }
   }

}
}
}
