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
#include <map>

// External includes
//

// Project includes
//

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
         IEnvironment();

         /**
          * @brief Destructor
          */
         virtual ~IEnvironment();

         /**
          * @brief Get the local ID
          */
         virtual int id() const;

         /**
          * @brief Size of environement (i.e CPUs, threads, etc)
          */
         virtual int size() const;

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
         virtual void abort(const int code) = 0;

         /**
          * @brief Abort with message
          */
         virtual void abort(const std::string msg) = 0;

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
