/** 
 * @file MpiFramework.cpp
 * @brief Source of the implementation of an MPI framework
 */

// Configuration includes
//

// System includes
//
#include <mpi.h>

// External includes
//

// Class include
//
#include "QuICC/Framework/MpiFramework.hpp"

// Project includes
//
#include "QuICC/Precision.hpp"
#include "QuICC/QuICCEnv.hpp"

namespace QuICC {

   void MpiFramework::initTransformComm(const int size)
   {
      MpiFramework::mTransformCpus.reserve(size);
      MpiFramework::mTransformComms.reserve(size);
      MpiFramework::mTransformIds.reserve(size);
   }

   void MpiFramework::addTransformComm(const ArrayI& ids)
   {
      // Store ranks of CPUs in group
      MpiFramework::mTransformCpus.push_back(ids);

      // MPI error code
      int ierr;

      // Split world communicator into sub groups
      MpiFramework::mTransformComms.push_back(MPI_Comm());
      int groupId = ids(0);
      ierr = MPI_Comm_split(MPI_COMM_WORLD, groupId, QuICCEnv().id(), &MpiFramework::mTransformComms.back());
      QuICCEnv().check(ierr, 913);

      // Make sure all communicators are real
      if(MpiFramework::mTransformComms.back() == MPI_COMM_NULL)
      {
         QuICCEnv().abort("MPI sub-communicators were not setup properly");
      }

      // Get rank in sub group
      int rank;
      ierr = MPI_Comm_rank(MpiFramework::mTransformComms.back(), &rank);
      QuICCEnv().check(ierr, 915);
      MpiFramework::mTransformIds.push_back(rank);

      // Check newly created communicator
      MpiFramework::checkTransformComm(MpiFramework::mTransformComms.size()-1, 111);
   }

   void MpiFramework::checkTransformComm(const int traId, const int code)
   {
      for(int i = 0; i < MpiFramework::transformCpus(traId).size(); ++i)
      {
         int rank = -1;
         if(MpiFramework::transformCpus(traId)(i) == QuICCEnv().id())
         {
            rank = QuICCEnv().id();

            MPI_Bcast(&rank, 1, MPI_INT, i, MpiFramework::transformComm(traId));

         } else
         {
            MPI_Bcast(&rank, 1, MPI_INT, i, MpiFramework::transformComm(traId));
         }

         if(rank != MpiFramework::transformCpus(traId)(i))
         {
            QuICCEnv().abort(code);
         }
      }
   }

   void MpiFramework::syncTransform(const int traId)
   {
      MPI_Barrier(MpiFramework::mTransformComms.at(traId));
   }

   int MpiFramework::transformId(const int traId)
   {
      return MpiFramework::mTransformIds.at(traId);
   }

   const ArrayI& MpiFramework::transformCpus(const int traId)
   {
      return MpiFramework::mTransformCpus.at(traId);
   }

   MPI_Comm MpiFramework::transformComm(const int traId)
   {
      return MpiFramework::mTransformComms.at(traId);
   }

   void MpiFramework::syncSubComm(const MpiFramework::SubCommId id, const int idx)
   {
      assert(MpiFramework::mSubComm.find(id) != MpiFramework::mSubComm.end());
      assert(MpiFramework::mSubComm.find(id)->second.size() > static_cast<size_t>(idx));

      // Create MPI barrier to force synchronisation
      MPI_Barrier(MpiFramework::mSubComm.find(id)->second.at(idx));
   }

   MPI_Group MpiFramework::getSubGroup(const MpiFramework::SubCommId id, const int idx)
   {
      assert(MpiFramework::mSubGroup.find(id) != MpiFramework::mSubGroup.end());

      return MpiFramework::mSubGroup.find(id)->second.at(idx);
   }

   MPI_Comm MpiFramework::getSubComm(const MpiFramework::SubCommId id, const int idx)
   {
      assert(MpiFramework::mSubComm.find(id) != MpiFramework::mSubComm.end());

      return MpiFramework::mSubComm.find(id)->second.at(idx);
   }

   void MpiFramework::finalize()
   {
      // Make sure all finished and are synchronised
      MPI_Barrier(MPI_COMM_WORLD);

      // Free communicators
      for(unsigned int i = 0; i < MpiFramework::mTransformComms.size(); i++)
      {
         MPI_Comm_free(&MpiFramework::mTransformComms.at(i));
      }

      // Free sub communicators
      for(auto subCommIt = MpiFramework::mSubComm.begin(); subCommIt != MpiFramework::mSubComm.end(); ++subCommIt)
      {
         for(unsigned int i = 0; i < subCommIt->second.size(); i++)
         {
            MPI_Comm_free(&subCommIt->second.at(i));
         }
      }

      // Free sub groups
      for(auto subGroupIt = MpiFramework::mSubGroup.begin(); subGroupIt != MpiFramework::mSubGroup.end(); ++subGroupIt)
      {
         for(unsigned int i = 0; i < subGroupIt->second.size(); i++)
         {
            MPI_Group_free(&subGroupIt->second.at(i));
         }
      }

      // Make sure all finished and are synchronised
      MPI_Barrier(MPI_COMM_WORLD);
   }

   void MpiFramework::initSubComm(const MpiFramework::SubCommId id, const int size)
   {
      // Add group and communicator
      mSubGroup.insert(std::make_pair(id, std::vector<MPI_Group>()));
      mSubComm.insert(std::make_pair(id, std::vector<MPI_Comm>()));

      for(int i = 0; i < size; i++)
      {
         mSubGroup.find(id)->second.push_back(MPI_GROUP_NULL);
         mSubComm.find(id)->second.push_back(MPI_COMM_NULL);
      }
   }

   void MpiFramework::setSubComm(const MpiFramework::SubCommId id, const int idx, const std::set<int>& ranks)
   {
      assert(ranks.size() > 0);
      assert(mSubGroup.find(id) != mSubGroup.end());
      assert(mSubComm.find(id) != mSubComm.end());

      // MPI error code
      int ierr;

      // Get world group
      MPI_Group   world;
      ierr = MPI_Comm_group(MPI_COMM_WORLD, &world);
      QuICCEnv().check(ierr, 921);

      MPI_Group group;
      MPI_Comm comm;

      // Create array of ranks
      ArrayI curRanks(ranks.size());
      int j = 0;
      for(auto sIt = ranks.begin(); sIt != ranks.end(); ++sIt)
      {
         curRanks(j) = *sIt;
         j++;
      }

      // Create sub group
      ierr = MPI_Group_incl(world, curRanks.size(), curRanks.data(), &group);
      QuICCEnv().check(ierr, 922);
      // Create sub communicator
      ierr = MPI_Comm_create(MPI_COMM_WORLD, group, &comm);
      QuICCEnv().check(ierr, 923);

      if(comm != MPI_COMM_NULL)
      {
         assert(mSubGroup.find(id)->second.size() > static_cast<size_t>(idx));
         assert(mSubComm.find(id)->second.size() > static_cast<size_t>(idx));

         mSubGroup.find(SPECTRAL)->second.at(idx) = group;
         mSubComm.find(SPECTRAL)->second.at(idx) = comm;
      }

      QuICCEnv().synchronize();
   }

   void MpiFramework::sleep(const MHDFloat seconds)
   {
      // Sleep MPI process for given amount of seconds
      MHDFloat start = MPI_Wtime();
      MHDFloat tmp = 0.0;
      while(MPI_Wtime() - start < seconds)
      {
         // Do Nothing ... just waiting
         tmp += 1.0;
      }
   }

   std::vector<int>  MpiFramework::mTransformIds = std::vector<int>();

   std::vector<ArrayI>  MpiFramework::mTransformCpus = std::vector<ArrayI>();

   std::vector<MPI_Comm>  MpiFramework::mTransformComms = std::vector<MPI_Comm>();

   std::map<MpiFramework::SubCommId,std::vector<MPI_Group> > MpiFramework::mSubGroup = std::map<MpiFramework::SubCommId,std::vector<MPI_Group> >();

   std::map<MpiFramework::SubCommId,std::vector<MPI_Comm> > MpiFramework::mSubComm = std::map<MpiFramework::SubCommId,std::vector<MPI_Comm> >();

}
