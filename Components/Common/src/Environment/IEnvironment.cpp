/** 
 * @file IEnvironment.cpp
 * @brief Source of generic environment
 */

// Configuration includes
//
#include "QuICC/Debug/DebuggerMacro.h"

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Environment/IEnvironment.hpp"

// Project includes
//

namespace QuICC {

namespace Environment {

   int IEnvironment::mIoRank = -99;

   int IEnvironment::mSize = -99;

   int IEnvironment::mId = -99;

   IEnvironment::IEnvironment()
   {
   }

   IEnvironment::~IEnvironment()
   {
   }

   int IEnvironment::ioRank() const
   {
      return this->mIoRank;
   }

   bool IEnvironment::allowsIO() const
   {
      return (this->ioRank() == this->id());
   }

   void IEnvironment::checkEnvironment(const int size)
   {
      // Check that the size was set
      if(this->size() < 1)
      {
         this->abort("Environment contains has negative size!");
      }

      // Check that the local ID was set
      if(this->id() < 0)
      {
         this->abort("Environment local ID was not set!");
      }

      // Check that IO rank was set
      if(this->ioRank() < 0)
      {
         this->abort("Environment IO rank was not set!");
      }

      // Check compatibility between requested cores and setup cores
      if(size != this->size())
      {
         this->abort("Environment parameters and setup have conflicting sizes: "
            +std::to_string(size)+" vs "+std::to_string(this->size()));
      }

      // Check profiler registration
      this->checkRegister<Register::Profiler>("Profiler");

      // Check multiple precision registration
      this->checkRegister<Register::Precision>("Multiple precision");

      // Check HDF5 registration
      this->checkRegister<Register::Hdf5>("HDF5");

      // Check MPI registration
      this->checkRegister<Register::Mpi>("MPI");
   }

   void IEnvironment::init()
   {
      DebuggerMacro_msg("Running environment initializers", 0);

      for(auto it = IEnvironment::initializers().begin(); it != IEnvironment::initializers().end(); ++it)
      {
         it->second();
      }
   }

   void IEnvironment::finalize()
   {
      DebuggerMacro_msg("Running environment finalizers", 0);

      for(auto it = IEnvironment::finalizers().rbegin(); it != IEnvironment::finalizers().rend(); ++it)
      {
         it->second();
      }
   }

   int IEnvironment::registerInitializer(const int id, std::function<void()> func)
   {
      if(IEnvironment::initializers().count(id) != 0)
      {
         throw std::logic_error("ID for Environment initialization already in use! Check priorities");
      }

      IEnvironment::initializers().insert(std::make_pair(id, func));

      return id;
   }

   int IEnvironment::registerFinalizer(const int id, std::function<void()> func)
   {
      if(IEnvironment::finalizers().count(id) != 0)
      {
         throw std::logic_error("ID for Environment finalization already in use! Check priorities");
      }

      IEnvironment::finalizers().insert(std::make_pair(id, func));

      return id;
   }

   std::map<int,std::function<void()> >& IEnvironment::initializers()
   {
      static std::map<int,std::function<void()> > inits;

      return inits;
   }

   std::map<int,std::function<void()> >& IEnvironment::finalizers()
   {
      static std::map<int,std::function<void()> > finals;

      return finals;
   }

}
}
