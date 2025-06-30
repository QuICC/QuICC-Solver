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
      : useDefault(true), dumpData(false), timeOnly(false), type(TestType::PROJECTOR), ulp(11), np(0), rank(0), iter(1)
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

   void TestArgs::clear()
   {
      this->useDefault = true;
      this->dumpData = false;
      this->timeOnly = false;
      this->type = TestType::PROJECTOR;
      this->ulp = 11;
      this->np = 0;
      this->rank = 0;
      this->iter = 1;
      this->params.clear();
   }

}
}
}
