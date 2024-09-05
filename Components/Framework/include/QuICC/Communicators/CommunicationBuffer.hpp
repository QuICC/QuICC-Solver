/**
 * @file CommunicationBuffer.hpp
 * @brief Implementation of a "raw" communication buffer
 */

#ifndef QUICC_PARALLEL_COMMUNICATIONBUFFER_HPP
#define QUICC_PARALLEL_COMMUNICATIONBUFFER_HPP

// System includes
//
#include <cassert>
#include <vector>

// Project includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"
#include "Types/BasicTypes.hpp"

namespace QuICC {

namespace Parallel {

/**
 *
 * @brief Implementation of a "raw" communication buffer
 */
template <typename TData> class CommunicationBuffer
{
public:
   /// Typedef for buffer datatype
   typedef TData DataType;

   /**
    * @brief Constructor
    */
   CommunicationBuffer();

   /**
    * @brief Destructor
    */
   ~CommunicationBuffer();

   /**
    * @brief Allocate buffers
    *
    * @param sizes      Sizes per node
    * @param maxPacks   Maximum number of packs
    */
   void allocate(const std::vector<int>& sizes, const int maxPacks);

   /**
    * @brief Allocate buffers to fit both requested sizes
    *
    * @param aSizes      Sizes per node
    * @param maxAPacks   Maximum number of packs
    * @param bSizes      Sizes per node
    * @param maxBPacks   Maximum number of packs
    */
   void allocateMax(const std::vector<int>& aSizes, const int maxAPacks,
      const std::vector<int>& bSizes, const int maxBPacks);

   /**
    * @brief Get pointer to raw buffer storage
    */
   TData* data();

   /**
    * @brief Get pointer to raw sub buffer start
    */
   TData* at(const int id);

   /**
    * @brief Get zero positions for sub buffer
    */
   int* zero();

   /**
    * @brief Get current position in sub buffer
    */
   int& pos(const int id);

   /**
    * @brief Reset positions in sub buffers
    */
   void resetPositions();

   /**
    * @brief Get total available storage
    */
   int total() const;

   /**
    * @brief Get the memory requirements
    */
   MHDFloat requiredStorage() const;

protected:
private:
   /**
    * @brief Cleanup the communication buffers
    */
   void cleanupBuffers();

   /**
    * @brief Total memory in buffer
    */
   int mTotal;

   /**
    * @brief MPI communication buffers
    */
   TData* mData;

   /**
    * @brief Start position for sub buffers
    */
   std::vector<int> mZero;

   /**
    * @brief Current position in sub buffers
    */
   std::vector<int> mPos;
};

template <typename TData> CommunicationBuffer<TData>::CommunicationBuffer() {}

template <typename TData> CommunicationBuffer<TData>::~CommunicationBuffer()
{
   // Cleanup the communication buffers
   this->cleanupBuffers();
}

template <typename TData>
void CommunicationBuffer<TData>::allocate(const std::vector<int>& sizes,
   const int maxPacks)
{
   // Create CPU group buffers
   this->mTotal = 0;
   for (auto it = sizes.cbegin(); it != sizes.cend(); ++it)
   {
      // Create zero position and initialize position
      this->mZero.push_back(this->mTotal);
      this->mPos.push_back(0);

      this->mTotal += (*it) * maxPacks;
   }

   // Allocate large buffer
   this->mData = new TData[this->mTotal];

#ifdef QUICC_STORAGEPROFILE
   MHDFloat mem = this->requiredStorage();

   StorageProfilerMacro_update(StorageProfilerMacro::MPI, mem);
   StorageProfilerMacro_update(StorageProfilerMacro::MPIBUFFERS, mem);
#endif // QUICC_STORAGEPROFILE
}

template <typename TData>
void CommunicationBuffer<TData>::allocateMax(const std::vector<int>& aSizes,
   const int maxAPacks, const std::vector<int>& bSizes, const int maxBPacks)
{
   // Create CPU group buffers
   this->mTotal = 0;
   for (std::size_t id = 0; id < std::min(aSizes.size(), bSizes.size()); ++id)
   {
      // Create zero position and initialize position
      this->mZero.push_back(this->mTotal);
      this->mPos.push_back(0);

      this->mTotal +=
         std::max(maxAPacks * aSizes.at(id), maxBPacks * bSizes.at(id));
   }

   // Deal with different number of CPUs in groups
   if (aSizes.size() > bSizes.size())
   {
      for (std::size_t id = bSizes.size(); id < aSizes.size(); ++id)
      {
         // Create zero position and initialize position
         this->mZero.push_back(this->mTotal);
         this->mPos.push_back(0);

         this->mTotal += maxAPacks * aSizes.at(id);
      }
   }
   else if (bSizes.size() > aSizes.size())
   {
      for (std::size_t id = aSizes.size(); id < bSizes.size(); ++id)
      {
         // Create zero position and initialize position
         this->mZero.push_back(this->mTotal);
         this->mPos.push_back(0);

         this->mTotal += maxBPacks * bSizes.at(id);
      }
   }

   this->mData = new TData[this->mTotal];

#ifdef QUICC_STORAGEPROFILE
   MHDFloat mem = this->requiredStorage();

   StorageProfilerMacro_update(StorageProfilerMacro::MPI, mem);
   StorageProfilerMacro_update(StorageProfilerMacro::MPIBUFFERS, mem);
#endif // QUICC_STORAGEPROFILE
}

template <typename TData> inline int CommunicationBuffer<TData>::total() const
{
   return this->mTotal;
}

template <typename TData> inline TData* CommunicationBuffer<TData>::data()
{
   assert(this->mData != nullptr);

   return this->mData;
}

template <typename TData> inline int* CommunicationBuffer<TData>::zero()
{
   assert(this->mZero.size() > 0);

   return &this->mZero[0];
}

template <typename TData>
inline TData* CommunicationBuffer<TData>::at(const int id)
{
   assert(this->mData != nullptr);

   return (this->mData + this->mZero.at(id));
}

template <typename TData> void CommunicationBuffer<TData>::resetPositions()
{
   for (auto it = this->mPos.begin(); it != this->mPos.end(); ++it)
   {
      *it = 0;
   }
}

template <typename TData>
inline int& CommunicationBuffer<TData>::pos(const int id)
{
   // Safety assert
   assert(this->mPos.size() > static_cast<std::size_t>(id));

   return this->mPos.at(id);
}

template <typename TData> void CommunicationBuffer<TData>::cleanupBuffers()
{
   // Free the buffers memory
   delete[] this->mData;
}

template <typename TData>
MHDFloat CommunicationBuffer<TData>::requiredStorage() const
{
   MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
   mem += Debug::MemorySize<TData>::BYTES * this->mTotal;
#endif // QUICC_STORAGEPROFILE

   return mem;
}

} // namespace Parallel
} // namespace QuICC

#endif // QUICC_PARALLEL_COMMUNICATIONBUFFER_HPP
