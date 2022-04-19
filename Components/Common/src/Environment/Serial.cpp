/** 
 * @file Serial.cpp
 * @brief Source of implementation of serial environment
 */

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Class include
//
#include "QuICC/Environment/Serial.hpp"

// Project includes
//
#include "QuICC/Precision.hpp"

namespace QuICC {

namespace Environment {

   Serial::Serial()
   {
   }

   Serial::~Serial()
   {
   }

   void Serial::init()
   {
      IEnvironment::init();

      // Set ID
      Serial::mId = 0;

      // Set IO rank
      Serial::mIoRank = Serial::mId;
   }

   void Serial::setup(const int size)
   {
      // Set the number of CPUs
      Serial::mSize = size;

      // Check that the environment was setup correctly
      this->checkEnvironment(1);
   }

   void Serial::synchronize()
   {
      // Nothing to be done in serial environment
   }

   void Serial::check(const int ierr, const int code)
   {
      if(ierr != 0)
      {
         this->abort(code);
      }
   }

   void Serial::abort(const int code)
   {
      throw std::logic_error("Aborted with error code: " + std::to_string(code));
   }

   void Serial::abort(const std::string msg)
   {
      throw std::logic_error("Aborted: " + msg);
   }

   void Serial::finalize()
   {
      IEnvironment::finalize();

   }

}
}
