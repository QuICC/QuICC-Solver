/**
 * @file Mpi.hpp
 * @brief Implementation of the MPI environment
 */

#ifdef QUICC_MPI
#ifndef QUICC_ENVIRONMENT_MPI_HPP
#define QUICC_ENVIRONMENT_MPI_HPP

// System includes
//
#include <mpi.h>

// External includes
//

// Project includes
//
#include "QuICC/Environment/IEnvironment.hpp"

namespace QuICC {

namespace Environment {

   /**
    * @brief This class defines the framework for an MPI Parallel code
    */
   class Mpi: public IEnvironment
   {
      public:
         /**
          * @brief Constructor
          */
         Mpi();

         /**
          * @brief Destructor
          */
         ~Mpi();

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
          * @brief Check error code for success
          */
         virtual void check(const int ierr, const std::string code) override;

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
          * @brief setup gdb hook
          */
         void gdbHook();

         /**
          * @brief size of the communicator
          */
         static int mCommSize;
   };

}
}

#endif // QUICC_ENVIRONMENT_MPI_HPP
#endif // QUICC_MPI
