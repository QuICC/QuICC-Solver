/**
 * @file MpiConverter.hpp
 * @brief Implementation of the MPI data converter
 */

#ifndef QUICC_PARALLEL_MPICONVERTER_HPP
#define QUICC_PARALLEL_MPICONVERTER_HPP

// Configuration includes
//
#include "QuICC/Debug/DebuggerMacro.h"

// System includes
//
#include <memory>
#include <type_traits>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Communicators/Converters/MpiConverterImpl.hpp"
#include "QuICC/Communicators/Converters/MpiConverterTools.hpp"
#include "QuICC/StorageProviders/DynamicPairProvider.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "QuICC/Communicators/Converters/IIndexConv.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of the MPI data converter.
    */
   class MpiConverter: public Framework::Selector::MpiConverterImpl
   {
      public:
         /// Typedef for forward datatype
         typedef Framework::Selector::MpiConverterImpl::RealFwdData RealFwdData;

         /// Typedef for forward datatype
         typedef Framework::Selector::MpiConverterImpl::ComplexFwdData ComplexFwdData;

         /// Typedef for backward datatype
         typedef Framework::Selector::MpiConverterImpl::RealBwdData RealBwdData;

         /// Typedef for backward datatype
         typedef Framework::Selector::MpiConverterImpl::ComplexBwdData ComplexBwdData;

         /**
          * @brief Constructor
          */
         MpiConverter();

         /**
          * @brief Destructor
          */
         virtual ~MpiConverter();

         /**
          * @brief Initialise the packs
          *
          * @param spRes      Shared Resolution
          * @param fwdDim     Dimension index for forward transform
          * @param fwdTmp     TFwd temporary
          * @param bwdTmp     TBwd temporary
          * @param fwdPacks   Array of possible pack sizes for forward transform
          * @param bwdPacks   Array of possible pack sizes for backward transform
          */
         template <typename TFwd,typename TBwd> void init(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TFwd& fwdTmp, TBwd& bwdTmp, const ArrayI& fwdPacks, const ArrayI& bwdPacks, std::shared_ptr<IIndexConv> spIdxConv);

         /**
          * @brief Convert real data from Fwd to Bwd
          */
         virtual void convertFwd(const RealFwdData &in, DynamicPairProvider &storage) override;

         /**
          * @brief Convert complex data from Fwd to Bwd
          */
         virtual void convertFwd(const ComplexFwdData &in, DynamicPairProvider &storage) override;

         /**
          * @brief Convert real data from Bwd to Fwd
          */
         virtual void convertBwd(const RealBwdData &in, DynamicPairProvider &storage) override;

         /**
          * @brief Convert complex data from Bwd to Fwd
          */
         virtual void convertBwd(const ComplexBwdData &in, DynamicPairProvider &storage) override;

         /**
          * @brief Get the converted real data from Bwd to Fwd
          */
         virtual void getFwd(RealFwdData *& pOut, DynamicPairProvider &storage) override;

         /**
          * @brief Get the converted complex data from Bwd to Fwd
          */
         virtual void getFwd(ComplexFwdData *& pOut, DynamicPairProvider &storage) override;

         /**
          * @brief Get the converted real data from Fwd to Bwd
          */
         virtual void getBwd(RealBwdData *& pOut, DynamicPairProvider &storage) override;

         /**
          * @brief Get the converted complex data from Fwd to Bwd
          */
         virtual void getBwd(ComplexBwdData *& pOut, DynamicPairProvider &storage) override;

         /**
         * @brief Do storage profiling
         */
         virtual void profileStorage() const override;

      protected:

      private:
         /**
          * @brief Convert real data from Fwd to Bwd
          */
         template <typename T> void convertFwdImpl(const T& in, DynamicPairProvider &storage);

         /**
          * @brief Convert real data from Bwd to Fwd
          */
         template <typename T> void convertBwdImpl(const T& in, DynamicPairProvider &storage);

         /**
          * @brief Get the converted complex data from Bwd to Fwd
          */
         template <typename T> void getFwdImpl(T *& pOut, DynamicPairProvider &storage);

         /**
          * @brief Get the converted real data from Fwd to Bwd
          */
         template <typename T> void getBwdImpl(T *& pOut, DynamicPairProvider &storage);

         /**
          * @brief Initialise the datatypes for packing
          *
          * @param spRes   Shared Resolution
          * @param fTmp    TFwd temporary
          * @param bTmp    TBwd temporary
          */
         template <typename TFwd,typename TBwd> void initPackTypes(SharedResolution spRes, TFwd& fTmp, TBwd& bTmp);

         /**
          * @brief Initialise the sizes and CPU lists
          */
         void initLists();
   };

   template <typename T> void MpiConverter::convertFwdImpl(const T& in, DynamicPairProvider&)
   {
      // Send the data
      this->sendFwd(in);
   }

   template <typename T> void MpiConverter::convertBwdImpl(const T& in, DynamicPairProvider&)
   {
      // Send the data
      this->sendBwd(in);
   }

   template <typename T> void MpiConverter::getFwdImpl(T *& pData, DynamicPairProvider &storage)
   {
      if(pData == nullptr)
      {
         // Get storage for output value
         storage.provideFwd(pData);
      }

      // Receive converted data
      this->receiveFwd(*pData);
   }

   template <typename T> void MpiConverter::getBwdImpl(T *& pData, DynamicPairProvider &storage)
   {
      if(pData == nullptr)
      {
         // Get storage for output value
         storage.provideBwd(pData);
      }

      // Receive converted data
      this->receiveBwd(*pData);
   }

   template <typename TFwd,typename TBwd> void MpiConverter::init(SharedResolution spRes, const Dimensions::Transform::Id fwdDim, TFwd &fwdTmp, TBwd &bwdTmp, const ArrayI& fwdPacks, const ArrayI& bwdPacks, std::shared_ptr<IIndexConv> spIdxConv)
   {
      StageTimer stage;
      stage.start("Creating MPI packing datatypes",1);

      // Set dimensions
      this->mDimensions = spRes->sim().ss().dimension();

      // Store the possible pack sizes
      this->mForwardPacks = fwdPacks;
      this->mBackwardPacks = bwdPacks;

      // Store Transform ID
      this->mTraId = fwdDim;

      // Set index converter
      this->setIndexConverter(spIdxConv);

      // initialise the data types
      this->initPackTypes(spRes, fwdTmp, bwdTmp);

      stage.done();
      stage.start("Creating CPU and size lists",1);

      // initialise the size and CPU lists
      this->initLists();

      stage.done();
   }

   template <typename TFwd, typename TBwd> void MpiConverter::initPackTypes(SharedResolution spRes, TFwd &fTmp, TBwd &bTmp)
   {
      MpiConverterTools::CoordinateMap localFwdMap;
      MpiConverterTools::buildLocalFwdMap(localFwdMap, spRes, this->mTraId);
      MpiConverterTools::CoordinateMap localBwdMap;
      MpiConverterTools::buildLocalBwdMap(localBwdMap, spRes, this->mTraId, this->idxConv());

      // Loop over group cpus
      this->mBTypes.reserve(QuICCEnv().size(this->mTraId));
      this->mFTypes.reserve(QuICCEnv().size(this->mTraId));
      for(int id = 0; id < QuICCEnv().size(this->mTraId); id++)
      {
      	 // Synchronize
         QuICCEnv().synchronize(this->mTraId);

         // Create TBwd datatypes
         #if defined QUICC_MPIPACK_MANUAL
            this->mBTypes.push_back(MpiConverterTools::CoordinateVector());
            MpiConverterTools::buildDatatype<TBwd,Dimensions::Data::DATF1D>(this->mBTypes.back(), localBwdMap, spRes, this->mTraId, this->mTraId, bTmp, id);
         #else
            MPI_Datatype type = MpiConverterTools::buildDatatype<TBwd,Dimensions::Data::DATF1D>(localBwdMap, spRes, this->mTraId, this->mTraId, bTmp, id);
            this->mBTypes.push_back(type);
         #endif //defined QUICC_MPIPACK_MANUAL

      	 // Synchronize
         QuICCEnv().synchronize(this->mTraId);

         // Create TFwd datatypes
         #if defined QUICC_MPIPACK_MANUAL
            this->mFTypes.push_back(MpiConverterTools::CoordinateVector());
            MpiConverterTools::buildDatatype<TFwd,Dimensions::Data::DATB1D>(this->mFTypes.back(), localFwdMap, spRes, this->mTraId, Dimensions::jump(this->mTraId,1), fTmp, id);
         #else
            type = MpiConverterTools::buildDatatype<TFwd,Dimensions::Data::DATB1D>(localFwdMap, spRes, this->mTraId, Dimensions::jump(this->mTraId,1), fTmp, id);
            this->mFTypes.push_back(type);
         #endif //defined QUICC_MPIPACK_MANUAL
      }
   }

}
}

#endif // QUICC_PARALLEL_MPICONVERTER_HPP
