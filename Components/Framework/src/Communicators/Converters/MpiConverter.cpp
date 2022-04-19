/**
 * @file MpiConverter.cpp
 * @brief Source of the serial data converter
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Communicators/Converters/MpiConverter.hpp"

// Project includes
//

namespace QuICC {

namespace Parallel {

   MpiConverter::MpiConverter()
   {
   }

   MpiConverter::~MpiConverter()
   {
   }

   void MpiConverter::convertFwd(const RealFwdData& in, DynamicPairProvider& storage)
   {
      this->convertFwdImpl(in, storage);
   }

   void MpiConverter::convertFwd(const ComplexFwdData& in, DynamicPairProvider& storage)
   {
      this->convertFwdImpl(in, storage);
   }

   void MpiConverter::convertBwd(const RealBwdData& in, DynamicPairProvider& storage)
   {
      this->convertBwdImpl(in, storage);
   }

   void MpiConverter::convertBwd(const ComplexBwdData& in, DynamicPairProvider& storage)
   {
      this->convertBwdImpl(in, storage);
   }

   void MpiConverter::getFwd(RealFwdData *& pData, DynamicPairProvider  &storage)
   {
      this->getFwdImpl(pData, storage);
   }

   void MpiConverter::getFwd(ComplexFwdData *& pData, DynamicPairProvider  &storage)
   {
      this->getFwdImpl(pData, storage);
   }

   void MpiConverter::getBwd(RealBwdData *& pData, DynamicPairProvider  &storage)
   {
      this->getBwdImpl(pData, storage);
   }

   void MpiConverter::getBwd(ComplexBwdData *& pData, DynamicPairProvider  &storage)
   {
      this->getBwdImpl(pData, storage);
   }

   void MpiConverter::initLists()
   {
      int sze;
      std::vector<int> unusedF;
      std::vector<int> unusedB;

      // Loop over all CPUs
      for(int i = 0; i < MpiFramework::transformCpus(this->mTraId).size(); i++)
      {
         // Compute buffer sizes for F group
         #if defined QUICC_MPIPACK_MANUAL
            sze = this->mFTypes.at(i).size();
         #else
            MPI_Pack_size(1, this->mFTypes.at(i), MpiFramework::transformComm(this->mTraId), &sze);
         #endif //defined QUICC_MPIPACK_MANUAL
         if(sze != 0 || this->mNeedEmptyComm)
         {
            this->mFSizes.push_back(sze);
            this->mFCpuGroup.push_back(i);

         // Get a list of unused forward datatypes
         } else
         {
            unusedF.push_back(i);
         }

         // Compute buffer sizes for B group
         #if defined QUICC_MPIPACK_MANUAL
            sze = this->mBTypes.at(i).size();
         #else
            MPI_Pack_size(1, this->mBTypes.at(i), MpiFramework::transformComm(this->mTraId), &sze);
         #endif //defined QUICC_MPIPACK_MANUAL
         if(sze != 0 || this->mNeedEmptyComm)
         {
            this->mBSizes.push_back(sze);
            this->mBCpuGroup.push_back(i);

         // Get a list of unused backward datatypes
         } else
         {
            unusedB.push_back(i);
         }
      }

      // Erase the unused forward datatypes
      for(auto rit = unusedF.rbegin(); rit != unusedF.rend(); ++rit)
      {
         this->mFTypes.erase(this->mFTypes.begin() + *rit);
      }

      // Erase the unused backward datatypes
      for(auto rit = unusedB.rbegin(); rit != unusedB.rend(); ++rit)
      {
         this->mBTypes.erase(this->mBTypes.begin() + *rit);
      }
   }

   void MpiConverter::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      MHDFloat memComm = 0.0;

      // General communication storage
      memComm += 4.0*5.0;
      memComm += 4.0*(this->mForwardPacks.size() + this->mBackwardPacks.size());
      memComm += 4.0*(this->mFSizes.size() + this->nFCpu());
      memComm += 4.0*(this->mBSizes.size() + this->nBCpu());

      // Requests communication storage
      memComm += 0.0;

      StorageProfilerMacro_update(StorageProfilerMacro::MPICOMM, memComm);
      StorageProfilerMacro_update(StorageProfilerMacro::MPI, memComm);

      MHDFloat memTypes = 0.0;

      // MPI datatypes storage
      memTypes += 0.0;

      StorageProfilerMacro_update(StorageProfilerMacro::MPITYPES, memTypes);
      StorageProfilerMacro_update(StorageProfilerMacro::MPI, memTypes);
#endif // QUICC_STORAGEPROFILE
   }

}
}
