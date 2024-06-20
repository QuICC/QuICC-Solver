/**
 * @file IEnvironment.hpp
 * @brief Generic interface for global coding environment stuff
 */

#ifndef QUICC_IENVIRONMENT_HPP
#define QUICC_IENVIRONMENT_HPP

// System includes
//
#include <cassert>
#include <string>
#include <vector>
#include <map>

// External includes
//

// Project includes
//
#include "Environment/Typedefs.hpp"

namespace QuICC {

namespace Environment {

   /**
    * @brief Generic interface for global coding environment stuff
    */
   class IEnvironment
   {
      public:
         /**
          * @brief Constructor
          */
         IEnvironment() = default;

         /**
          * @brief Destructor
          */
         virtual ~IEnvironment() = default;

         /**
          * @brief Get the local ID
          */
         virtual int id() const;

         /**
          * @brief Get rank in sub-communicator
          */
         virtual int id(const std::size_t id) const = 0;

         /**
          * @brief Size of environement (i.e CPUs, threads, etc)
          */
         virtual int size() const;

         /**
          * @brief Size of environement of sub-communicator (i.e CPUs, threads, etc)
          */
         virtual int size(const std::size_t id) const = 0;

         /**
          * @brief Get the IO rank
          */
         virtual int ioRank() const;

         /**
          * @brief Check if local id is allowed to do IO
          */
         virtual bool allowsIO() const;

         /**
          * @brief Setup the environment
          */
         virtual void setup(const int size) = 0;

         /**
          * @brief Synchronise
          */
         virtual void synchronize() = 0;

         /**
          * @brief Synchronize sub-communicator
          */
         virtual void synchronize(const std::size_t id) const = 0;

         /**
          * @brief Check error code for success
          */
         virtual void check(const int ierr, const int code) = 0;

         /**
          * @brief Check error code for success
          */
         virtual void check(const int ierr, const std::string msg) = 0;

         /**
          * @brief Abort with error code
          */
         virtual void abort(const int code) const = 0;

         /**
          * @brief Abort with message
          */
         virtual void abort(const std::string msg) const = 0;

         /**
          * @brief Add CPU group IDs
          */
         virtual void addCommunicator(const std::size_t id, const std::vector<int>& ids) = 0;

         /**
          * @brief Get WORLD ranks in sub-communicator group
          */
         virtual const std::vector<int>& groupIds(const std::size_t id) const = 0;

         /**
          * @brief Get MPI sub-communicator
          */
         virtual CommType comm(const std::size_t id) const = 0;

      protected:

         /**
          * @brief Number of computation cores
          */
         static int mSize;

         /**
          * @brief ID of the current node (ie MPI rank)
          */
         static int mId;

         /**
          * @brief ID allowed to do basic serial IO
          */
         static int mIoRank;

         /**
          * @brief Check setup of environment
          *
          * @param size Configured size
          */
         virtual void checkEnvironment(const int size);
   };

   inline int IEnvironment::size() const
   {
      // Safety assert to avoid uninitialised use
      assert(IEnvironment::mSize > 0);

      return IEnvironment::mSize;
   }

   inline int IEnvironment::id() const
   {
      // Safety assert to avoid uninitialised use
      assert(IEnvironment::mId >= 0);

      return IEnvironment::mId;
   }

}
}

#endif // QUICC_IENVIRONMENT_HPP
