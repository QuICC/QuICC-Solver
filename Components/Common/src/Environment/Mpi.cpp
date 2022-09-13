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
#include <climits>
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
      MPI_Init(nullptr, nullptr);

      // Set MPI rank of local CPU
      MPI_Comm_rank(MPI_COMM_WORLD, &Mpi::mId);

      // Get communnicator size
      MPI_Comm_size(MPI_COMM_WORLD, &Mpi::mCommSize);

      // Set IO rank
      Mpi::mIoRank = 0;

      // Set gdb hook
      gdbHook();
   }

   Mpi::~Mpi()
   {
      MPI_Finalize();
   }

   void Mpi::setup(const int size)
   {
      // Set the number of ranks
      Mpi::mSize = size;

      // Check that the environment was setup correctly
      this->checkEnvironment(mCommSize);
   }

   void Mpi::synchronize()
   {
      MPI_Barrier(MPI_COMM_WORLD);
   }

   void Mpi::check(const int ierr, const int code)
   {
      #ifndef NDEBUG
         if(ierr != MPI_SUCCESS)
         {
            this->abort(code);
         }
      #endif
   }

   void Mpi::check(const int ierr, const std::string msg)
   {
      #ifndef NDEBUG
         if(ierr != MPI_SUCCESS)
         {
            this->abort(msg);
         }
      #endif
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
      for(int i = 0; i < 100; ++i)
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
         std::cerr << "Aborted: " + msg << std::endl;
         using namespace std::chrono_literals;
         std::this_thread::sleep_for(1ms);
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

   void Mpi::gdbHook()
   {
      const char* env = std::getenv("QUICC_GDB_HOOK");
      if (env)
      {
         // Try to get rank from env
         int gdbRank;
         try
         {
            gdbRank = std::stoi(env) % mCommSize;
         }
         catch (...)
         {
            this->abort("gdb rank needs to be a non-negative integer");
         }

         if(gdbRank < 0)
         {
            this->abort("gdb rank needs to be non-negative");
         }

         // gdb rank wait for gdb connection
         if(Mpi::mId == gdbRank)
         {
            using namespace std::chrono_literals;
            volatile int i = 0;
            char hostname[HOST_NAME_MAX];
            gethostname(hostname, HOST_NAME_MAX);
            std::cerr << "PID " << getpid() << " on " << hostname
               << " with rank " << gdbRank
               << " ready to attach" << std::endl;
            while (0 == i)
            std::this_thread::sleep_for(1s);
         }
      }
      // Wait for everyone
      this->synchronize();
   }

   // Static init
   int Mpi::mCommSize = -99;

} // namespace Enviroment
} // namespace QuICC
