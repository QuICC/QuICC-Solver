/**
 * @file Mpi.cpp
 * @brief Source of the implementation of an MPI environment
 */

// Configuration includes
//

// System includes
//
#include <chrono>
#include <climits>
#include <iostream>
#include <regex>
#include <thread>
extern "C" {
#include <unistd.h>
}

// External includes
//

// Class include
//
#include "Environment/Mpi.hpp"

namespace QuICC {

namespace Environment {

Mpi::Mpi()
{
   MPI_Init(nullptr, nullptr);

   // Set MPI rank of local CPU
   MPI_Comm_rank(MPI_COMM_WORLD, &Mpi::mId);

   // Get communicator size
   MPI_Comm_size(MPI_COMM_WORLD, &Mpi::mCommSize);

   // Set IO rank
   Mpi::mIoRank = 0;

   // Set gdb hook
   gdbHook();
}

Mpi::~Mpi()
{
   this->synchronize();

   // Free communicators
   for (auto&& c: this->mComms)
   {
      MPI_Comm_free(&c.second);
   }

   this->synchronize();

   MPI_Finalize();
}

void Mpi::setup(const int size)
{
   // Set the number of ranks
   Mpi::mSize = size;

   // Check that the environment was setup correctly
   this->checkEnvironment(mCommSize);
}

int Mpi::id(const std::size_t id) const
{
   assert(this->mGroupId.count(id) > 0);

   return this->mGroupId.at(id);
}

int Mpi::size(const std::size_t id) const
{
   assert(this->mGroupId.count(id) > 0);

   return this->mWorldIds.at(id).size();
}

void Mpi::synchronize()
{
   MPI_Barrier(MPI_COMM_WORLD);
}

void Mpi::synchronize(const std::size_t id) const
{
   assert(this->mComms.count(id) > 0);

   MPI_Barrier(this->mComms.at(id));
}

void Mpi::check(const int ierr, const int code)
{
#ifndef NDEBUG
   if (ierr != MPI_SUCCESS)
   {
      this->abort(code);
   }
#endif
}

void Mpi::check(const int ierr, const std::string msg)
{
#ifndef NDEBUG
   if (ierr != MPI_SUCCESS)
   {
      this->abort(msg);
   }
#endif
}

void Mpi::abort(const int code) const
{
   MPI_Abort(MPI_COMM_WORLD, code);
}

void Mpi::abort(const std::string msg) const
{
   MPI_Request request;
   MPI_Ibarrier(MPI_COMM_WORLD, &request);

   // let's wait a bit for everyone to join
   int request_complete;
   MPI_Status status;
   for (int i = 0; i < 100; ++i)
   {
      MPI_Test(&request, &request_complete, &status);
      if (request_complete)
      {
         break;
      }
      using namespace std::chrono_literals;
      std::this_thread::sleep_for(10ms);
   }

   // we tried...
   if (!request_complete)
   {
      std::cerr << "Aborted: " + msg << std::endl;
      using namespace std::chrono_literals;
      std::this_thread::sleep_for(1ms);
      MPI_Abort(MPI_COMM_WORLD, 666);
   }
   MPI_Wait(&request, MPI_STATUS_IGNORE);

   // everyone makes it here
   if (Mpi::mId == Mpi::mIoRank)
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

      if (gdbRank < 0)
      {
         this->abort("gdb rank needs to be non-negative");
      }

      // gdb rank wait for gdb connection
      if (Mpi::mId == gdbRank)
      {
         using namespace std::chrono_literals;
         volatile bool wait = true;
         char hostname[HOST_NAME_MAX];
         gethostname(hostname, HOST_NAME_MAX);
         std::cerr << "PID " << getpid() << " on " << hostname << " with rank "
                   << gdbRank << " ready to attach" << std::endl;
         while (wait == true)
         {
            std::this_thread::sleep_for(1s);
         }
      }
   }
   // Wait for everyone
   this->synchronize();
}

void Mpi::addCommunicator(const std::size_t id, const std::vector<int>& ids)
{
   // Store ranks of CPUs in group
   this->mWorldIds.emplace(id, ids);

   // MPI error code
   int ierr;

   // Split world communicator into sub groups
   this->mComms.emplace(id, MPI_Comm());
   int groupId = ids.at(0);
   ierr =
      MPI_Comm_split(MPI_COMM_WORLD, groupId, this->id(), &this->mComms.at(id));
   this->check(ierr, "Splitting of MPI_COMM_WORLD failed");

   // Make sure setup of communicators was sucessful
   if (this->mComms.at(id) == MPI_COMM_NULL)
   {
      this->abort("MPI sub-communicators were not setup properly");
   }

   // Get rank in sub group
   int rank;
   ierr = MPI_Comm_rank(this->mComms.at(id), &rank);
   this->check(ierr, "Getting local rank in sub-communicator failed");
   this->mGroupId.emplace(id, rank);

   // Check newly created communicator
   this->checkCommunicator(id);
}

void Mpi::checkCommunicator(const std::size_t id) const
{
   for (int i = 0; i < this->size(id); ++i)
   {
      int rank = -1;
      if (this->mWorldIds.at(id).at(i) == this->id())
      {
         rank = this->id();

         MPI_Bcast(&rank, 1, MPI_INT, i, this->mComms.at(id));
      }
      else
      {
         MPI_Bcast(&rank, 1, MPI_INT, i, this->mComms.at(id));
      }

      if (rank != this->mWorldIds.at(id).at(i))
      {
         this->abort("Checking of MPI sub-communicator ranks failed");
      }
   }
}

const std::vector<int>& Mpi::groupIds(const std::size_t id) const
{
   assert(this->mWorldIds.count(id) > 0);

   return this->mWorldIds.at(id);
}

MPI_Comm Mpi::comm(const std::size_t id) const
{
   assert(this->mComms.count(id) > 0);

   return this->mComms.at(id);
}

// Static init
int Mpi::mCommSize = -99;

} // namespace Environment
} // namespace QuICC
