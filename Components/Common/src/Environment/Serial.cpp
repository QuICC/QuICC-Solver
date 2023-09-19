/**
 * @file Serial.cpp
 * @brief Source of implementation of serial environment
 */

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Class include
//
#include "QuICC/Environment/Serial.hpp"

// Project includes
//
#include "Types/Precision.hpp"

namespace QuICC {

namespace Environment {

   Serial::Serial()
   {
      // Set ID
      Serial::mId = 0;

      // Set IO rank
      Serial::mIoRank = Serial::mId;
   }

   Serial::~Serial()
   {
   }

   int Serial::id(const std::size_t id) const
   {
      return 0;
   }

   int Serial::size(const std::size_t id) const
   {
      return 1;
   }

   void Serial::setup(const int size)
   {
      // Preconditions
      checkMpiLauncher();
      if(size != 1)
      {
         this->abort("Serial setup can only use 1 cpu!");
      }

      // Set the number of CPUs
      Serial::mSize = size;

      // Check that the environment was setup correctly
      this->checkEnvironment(1);
   }

   void Serial::synchronize()
   {
      // Nothing to be done in serial environment
   }

   void Serial::synchronize(const std::size_t id) const
   {
      // Nothing to be done in serial environment
   }

   void Serial::check(const int ierr, const int code)
   {
      #ifndef NDEBUG
         if(ierr != 0)
         {
            this->abort(code);
         }
      #endif
   }

   void Serial::check(const int ierr, const std::string msg)
   {
      #ifndef NDEBUG
         if(ierr != 0)
         {
            this->abort(msg);
         }
      #endif
   }

   void Serial::abort(const int code) const
   {
      throw std::logic_error("Aborted with error code: " + std::to_string(code));
   }

   void Serial::abort(const std::string msg) const
   {
      throw std::logic_error("Aborted: " + msg);
   }

   void Serial::addCommunicator(const std::size_t, const std::vector<int>&)
   {
      // Nothing to be done in serial environment
   }

   const std::vector<int>& Serial::groupIds(const std::size_t) const
   {
      static std::vector<int> empty;
      return empty;
   }

   int Serial::comm(const std::size_t id) const
   {
      return -1;
   }

   void Serial::checkMpiLauncher()
   {
      // Check Open MPI (guaranteed)
      auto opmiCheck = std::getenv("OMPI_COMM_WORLD_SIZE");
      if(opmiCheck != nullptr)
      {
         this->abort("Serial run launched with Open MPI!");
      }

      // Check MPICH (not guaranteed)
      auto mpichCheck1 = std::getenv("MPICH_RANK_REORDER_METHOD");
      auto mpichCheck2 = std::getenv("MPICH_MPIIO_HINTS");
      if(mpichCheck1 != nullptr || mpichCheck2 != nullptr)
      {
         this->abort("Serial run launched with MPICH!");
      }

      // Check srun
      auto srunCheck = std::getenv("SLURM_NTASKS");
      if(srunCheck != nullptr)
      {
         int nTasks = std::stoi(srunCheck);
         // srun and 1 task is needed for interactive jobs
         if (nTasks != 1)
         {
            this->abort("Serial run launched with srun!");
         }
      }

   }

} // namespace Environment
} // namespace QuICC
