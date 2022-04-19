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
#include <functional>

// External includes
//

// Project includes
//
#include "QuICC/Environment/Register/Profiler.hpp"
#include "QuICC/Environment/Register/Precision.hpp"
#include "QuICC/Environment/Register/Hdf5.hpp"
#include "QuICC/Environment/Register/Mpi.hpp"

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
          * @brief Initialise the environment
          */
         virtual void init();

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
          * @brief Abort with error code
          */
         virtual void abort(const int code) = 0;

         /**
          * @brief Abort with message
          */
         virtual void abort(const std::string msg) = 0;

         /**
          * @brief Finalise environment
          */
         virtual void finalize();

         /**
          * @brief Register an initializer
          */
         static int registerInitializer(const int id, std::function<void()> fct);

         /**
          * @brief Register an initializer
          */
         static int registerFinalizer(const int id, std::function<void()> fct);

         /**
          * @brief Map of initializers (library, components, etc)
          */
         static std::map<int,std::function<void()> >& initializers();

         /**
          * @brief Map of finalizers (library, components, etc)
          */
         static std::map<int,std::function<void()> >& finalizers();

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

      private:
         /**
          * @brief Check registered initializer/finalizer
          */
         template <typename T> void checkRegister(const std::string name);
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

   template <typename T> void IEnvironment::checkRegister(const std::string name)
   {
      if(T::isActive && IEnvironment::initializers().count(T::initId) != 1)
      {
         this->abort(name + " initializer was not registerd");
      }

      if(T::isActive && IEnvironment::finalizers().count(T::finalId) != 1)
      {
         this->abort(name + " finalizer was not registerd");
      }
   }

}
}

#endif // QUICC_IENVIRONMENT_HPP
