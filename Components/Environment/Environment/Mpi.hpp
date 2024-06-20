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
#include <vector>
#include <map>

// External includes
//

// Project includes
//
#include "Environment/IEnvironment.hpp"

namespace QuICC {

namespace Environment {

   /**
    * @brief This class defines the framework for an MPI Parallel code
    */
   class Mpi: public IEnvironment
   {
      public:
         /**
          * @name Enum for different sub communicators
          */
         enum SubCommId {
            /// Spectral space CPUs
            SPECTRAL,
         };

         /**
          * @brief Constructor
          */
         Mpi();

         /**
          * @brief Destructor
          */
         ~Mpi();

         using IEnvironment::id;

         /**
          * @brief Get rank in sub-communicator
          */
         int id(const std::size_t id) const final;

         using IEnvironment::size;

         /**
          * @brief Get size of sub-communicator
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
         void check(const int ierr, const std::string code) final;

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
          * @brief Check MPI sub-communicator
          */
         void checkCommunicator(const std::size_t id) const;

         /**
          * @brief setup gdb hook
          */
         void gdbHook();

         /**
          * @brief size of the communicator
          */
         static int mCommSize;

         /**
          * @brief Local CPU rank in transform group
          */
         std::map<std::size_t,int> mGroupId;

         /**
          * @brief IDs of the CPUs in transform communication groups
          */
         std::map<std::size_t,std::vector<int> > mWorldIds;

         /**
          * @brief MPI communicators of the CPUs in transform communication groups
          */
         std::map<std::size_t,MPI_Comm> mComms;
   };

}
}

#endif // QUICC_ENVIRONMENT_MPI_HPP
#endif // QUICC_MPI
