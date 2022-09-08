/**
 * @file Serial.hpp
 * @brief Implementation of the serial environement
 */

#ifndef QUICC_ENVIRONMENT_SERIAL_HPP
#define QUICC_ENVIRONMENT_SERIAL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Environment/IEnvironment.hpp"

namespace QuICC {

namespace Environment {

   /**
    * @brief Serial environment
    */
   class Serial: public IEnvironment
   {
      public:
         /**
          * @brief Constructor
          */
         Serial();

         /**
          * @brief Destructor
          */
         ~Serial();

         /**
          * @brief Setup the environment
          */
         virtual void setup(const int size) override;

         /**
          * @brief Synchronise
          */
         virtual void synchronize() override;

         /**
          * @brief Check error code for success
          */
         virtual void check(const int ierr, const int code) override;

         /**
          * @brief Abort with error code
          */
         virtual void abort(const int code) override;

         /**
          * @brief Abort with message
          */
         virtual void abort(const std::string msg) override;

      protected:

      private:

         /**
          * @brief Check for known mpi launchers
          */
         void checkMpiLauncher();
   };

}
}

#endif // QUICC_ENVIRONMENT_SERIAL_HPP
