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
#include "QuICC/TestSuite/Polynomial/TestArgs.hpp"

namespace QuICC {

namespace TestSuite {

namespace Polynomial {

   TestArgs::TestArgs()
      : useDefault(true), dumpData(false), type(TestType::MATRIX), specN(0), physN(0), ulp(11)
   {
   }

   void TestArgs::setType(const std::string& type)
   {
      if(type == "matrix")
      {
         this->type = TestType::MATRIX;
      }
      else if(type == "weighted")
      {
         this->type = TestType::WEIGHTED_MATRIX;
      }
      else if(type == "otf_inner")
      {
         this->type = TestType::OTF_INNER;
      }
      else if(type == "otf_outer")
      {
         this->type = TestType::OTF_OUTER;
      }
      else if(type == "otf_reduce")
      {
         this->type = TestType::OTF_REDUCE;
      }
      else if(type == "quadrature")
      {
         this->type = TestType::QUADRATURE;
      }
      else
      {
         throw std::logic_error("Unsupported test type!");
      }
   }

}
}
}
