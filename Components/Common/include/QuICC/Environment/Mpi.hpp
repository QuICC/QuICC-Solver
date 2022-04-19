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
          * @brief Initialise the environment
          */
         virtual void init() override;

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

         /**
          * @brief Finalise environment
          */
         virtual void finalize() override;

      protected:

      private:
   };

}
}

#endif // QUICC_ENVIRONMENT_MPI_HPP
#endif // QUICC_MPI
