/**
 * @file MpiConverterBase.hpp
 * @brief Templated implementation of the base of a MPI data converter
 */

#ifndef QUICC_PARALLEL_MPICONVERTERBASE_HPP
#define QUICC_PARALLEL_MPICONVERTERBASE_HPP

// Debug includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

// Configuration includes
//

// System includes
//
#include <cassert>
#include <set>
#include <map>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/MpiTypes.hpp"
#include "QuICC/Enums/TransformDirection.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Communicators/Converters/IConverter.hpp"
#include "QuICC/Communicators/Converters/MpiConverterTools.hpp"
#include "QuICC/Communicators/CommunicationBuffer.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Templated implementation of the base of a MPI data converter.
    */
   class MpiConverterBase: public IConverter
   {
      public:
         /// Typedef for forward datatype
         typedef IConverter::RealFwdData RealFwdData;

         /// Typedef for forward datatype
         typedef IConverter::ComplexFwdData ComplexFwdData;

         /// Typedef for backward datatype
         typedef IConverter::RealBwdData RealBwdData;

         /// Typedef for backward datatype
         typedef IConverter::ComplexBwdData ComplexBwdData;

         /// Typedef for communication buffer for forward data
         typedef CommunicationBuffer<char> FwdBufferType;

         /// Typedef for communication buffer for backward data
         typedef CommunicationBuffer<char> BwdBufferType;

         /// Typedef for shared communication buffer for forward data
         typedef std::shared_ptr<FwdBufferType> SharedFwdBufferType;

         /// Typedef for shared communication buffer for backward data
         typedef std::shared_ptr<BwdBufferType> SharedBwdBufferType;

         /**
          * @brief Constructor
          */
         MpiConverterBase();

         /**
          * @brief Destructor
          */
         virtual ~MpiConverterBase();

         /**
          * @brief Set the communication buffers
          *
          * @brief spFwd Forward communication buffers
          * @brief spBwd Backward communication buffers
          */
         void setBuffers(SharedFwdBufferType spFwd, SharedBwdBufferType spBwd);

         /**
          * @brief Get forward buffer sizes
          */
         const std::vector<int> & fwdSizes() const;

         /**
          * @brief Get backward buffer sizes
          */
         const std::vector<int> & bwdSizes() const;

         /**
         * @brief Do storage profiling
         */
         virtual void profileStorage() const override;

      protected:
         /**
          * @brief Reset Fwd buffer positions
          */
         void resetFwdPositions();

         /**
          * @brief Reset Bwd buffer positions
          */
         void resetBwdPositions();

         /**
          * @brief Size of the forward packet
          *
          * @param id ID of the node
          */
         int sizeFPacket(const int id) const;

         /**
          * @brief Size of the backward packet
          *
          * @param id ID of the node
          */
         int sizeBPacket(const int id) const;

         /**
          * @brief Get size of the forward CPU group
          */
         int nFCpu() const;

         /**
          * @brief Get size of the backward CPU group
          */
         int nBCpu() const;

         /**
          * @brief Get MPI rank of CPU from forward CPU group
          *
          * @param id CPU group id
          */
         int fCpu(const int id) const;

         /**
          * @brief Get MPI rank of CPU from backward CPU group
          *
          * @param id CPU group id
          */
         int bCpu(const int id) const;

         /**
          * @brief Keep empty communcations?
          */
         bool mNeedEmptyComm;

         /**
          * @brief Communication packs counter
          */
         int mPacks;

         /**
          * @brief Direction of operations
          */
         TransformDirection::Id   mDirection;

         /**
          * @brief Transform ID
          */
         Dimensions::Transform::Id mTraId;

         /**
          * @brief List of CPU ranks involved in the forward conversion
          */
         std::vector<int>  mFCpuGroup;

         /**
          * @brief List of CPU ranks involved in the backward conversion
          */
         std::vector<int>  mBCpuGroup;

         /**
          * @brief Forward communication
          */
         SharedFwdBufferType mspFBuffers;

         /**
          * @brief Backward communication buffer pointer
          */
         SharedBwdBufferType mspBBuffers;

         /**
          * @brief List of the forward buffer sizes
          */
         std::vector<int>  mFSizes;

         /**
          * @brief List of the backward buffer sizes
          */
         std::vector<int>  mBSizes;

         /**
          * @brief Possible forward transform packs
          */
         ArrayI   mForwardPacks;

         /**
          * @brief Possible backward transform packs
          */
         ArrayI   mBackwardPacks;

         #if defined QUICC_MPIPACK_MANUAL
            /**
             * @brief Storage for the forward datatypes
             */
            std::vector<MpiConverterTools::CoordinateVector>  mFTypes;

            /**
             * @brief Storage for the backward datatypes
             */
            std::vector<MpiConverterTools::CoordinateVector> mBTypes;

         #else
            /**
             * @brief Storage for the forward datatypes
             */
            std::vector<MPI_Datatype>  mFTypes;

            /**
             * @brief Storage for the backward datatypes
             */
            std::vector<MPI_Datatype> mBTypes;
         #endif //defined QUICC_MPIPACK_MANUAL

      private:
         /**
          * @brief Cleanup the data types
          */
         void cleanupTypes();
   };

   inline const std::vector<int>& MpiConverterBase::fwdSizes() const
   {
      return this->mFSizes;
   }

   inline const std::vector<int>& MpiConverterBase::bwdSizes() const
   {
      return this->mBSizes;
   }

   inline int MpiConverterBase::nFCpu() const
   {
      return this->mFCpuGroup.size();
   }

   inline int MpiConverterBase::nBCpu() const
   {
      return this->mBCpuGroup.size();
   }

   inline int MpiConverterBase::fCpu(const int id) const
   {
      return this->mFCpuGroup.at(id);
   }

   inline int MpiConverterBase::bCpu(const int id) const
   {
      return this->mBCpuGroup.at(id);
   }

   inline int MpiConverterBase::sizeFPacket(const int id) const
   {
      return this->mPacks*this->mFSizes.at(id);
   }

   inline int MpiConverterBase::sizeBPacket(const int id) const
   {
      return this->mPacks*this->mBSizes.at(id);
   }

}
}

#endif // QUICC_PARALLEL_MPICONVERTERBASE_HPP
