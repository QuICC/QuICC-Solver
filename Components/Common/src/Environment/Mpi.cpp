/** 
 * @file Mpi.cpp
 * @brief Source of the implementation of an MPI environment
 */

// Configuration includes
//

// System includes
//
#include <iostream>
#include <chrono>
#include <thread>
#include <regex>
extern "C"
{
   #include <unistd.h>
}

// External includes
//

// Class include
//
#include "QuICC/Environment/Mpi.hpp"

// Project includes
//

namespace QuICC {

namespace Environment {

   Mpi::Mpi()
   {
   }

   Mpi::~Mpi()
   {
   }

   void Mpi::init()
   {
      IEnvironment::init();

      // Set MPI rank of local CPU
      MPI_Comm_rank(MPI_COMM_WORLD, &Mpi::mId);

      // Set IO rank
      Mpi::mIoRank = 0;

      // Set gdb hook
      gdbHook();
   }

   void Mpi::setup(const int size)
   {
      // Set the number of ranks
      Mpi::mSize = size;

      // Get MPI size
      int nRank;
      MPI_Comm_size(MPI_COMM_WORLD, &nRank);

      // Check that the environment was setup correctly
      this->checkEnvironment(nRank);
   }

   void Mpi::synchronize()
   {
      MPI_Barrier(MPI_COMM_WORLD);
   }

   void Mpi::check(const int ierr, const int code)
   {
      if(ierr != MPI_SUCCESS)
      {
         this->abort(code);
      }
   }

   void Mpi::abort(const int code)
   {
      MPI_Abort(MPI_COMM_WORLD, code);
   }

   void Mpi::abort(const std::string msg)
   {
      MPI_Request request;
      MPI_Ibarrier(MPI_COMM_WORLD, &request);

      // let's wait a bit for everyone to join
      int request_complete;
      MPI_Status status;
      for(int i = 0; i < 10; ++i)
      {
         MPI_Test(&request, &request_complete, &status);
         if(request_complete)
         {
            break;
         }
         using namespace std::chrono_literals;
         std::this_thread::sleep_for(10ms);
      }

      // we tried...
      if(!request_complete)
      {
         MPI_Abort(MPI_COMM_WORLD, 666);
      }
      MPI_Wait(&request, MPI_STATUS_IGNORE);

      // everyone makes it here
      if(Mpi::mId == Mpi::mIoRank)
      {
         std::cerr << "Aborted: " + msg << std::endl;
      }
      // after all this effort let's make sure to print!
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
   }

   void Mpi::finalize()
   {
      // Make sure all finished and are synchronised
      MPI_Barrier(MPI_COMM_WORLD);

      IEnvironment::finalize();
   }

   void Mpi::gdbHook()
   {
      const char* env = std::getenv("QUICC_GDB_HOOK");
      if (env) {
         if(Mpi::mId == Mpi::mIoRank)
         {
            using namespace std::chrono_literals;
            volatile int i = 0;
            char hostname[HOST_NAME_MAX];
            gethostname(hostname, HOST_NAME_MAX);
            std::cerr << "PID " << getpid() << " on " << hostname
               << " ready to attach" << std::endl;
            while (0 == i)
            std::this_thread::sleep_for(1s);
         }
      }
   }

} // namespace Enviroment
} // namespace QuICC
