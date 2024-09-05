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
#include "Environment/IEnvironment.hpp"

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

         using IEnvironment::id;

         /**
          * @brief Get rank in sub-communicator
          */
         int id(const std::size_t id) const final;

         using IEnvironment::size;

         /**
          * @brief Size of environement of sub-communicator (i.e CPUs, threads, etc)
          */
         int size(const std::size_t id) const final;

         /**
          * @brief Setup the environment
          */
         void setup(const int size) final;

         /**
          * @brief Synchronise
          */
         void synchronize() final;

         /**
          * @brief Synchronize sub-communicator
          */
         void synchronize(const std::size_t id) const final;

         /**
          * @brief Check error code for success
          */
         void check(const int ierr, const int code) final;

         /**
          * @brief Check error code for success
          */
         void check(const int ierr, const std::string msg) final;

         /**
          * @brief Abort with error code
          */
         void abort(const int code) const final;

         /**
          * @brief Abort with message
          */
         void abort(const std::string msg) const final;

         /**
          * @brief Add CPU group IDs
          */
         void addCommunicator(const std::size_t id, const std::vector<int>& ids) final;

         /**
          * @brief Get WORLD ranks in sub-communicator group
          */
         const std::vector<int>& groupIds(const std::size_t id) const final;

         /**
          * @brief Get MPI sub-communicator
          */
         CommType comm(const std::size_t id) const final;

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
