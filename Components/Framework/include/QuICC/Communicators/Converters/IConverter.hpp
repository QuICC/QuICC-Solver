/**
 * @file IConverter.hpp
 * @brief Implementation of the interface for a data converter
 */

#ifndef QUICC_PARALLEL_ICONVERTER_HPP
#define QUICC_PARALLEL_ICONVERTER_HPP

// Configuration includes
//

// System includes
//
#include <memory>
#include <type_traits>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/TransformDirection.hpp"
#include "QuICC/StorageProviders/DynamicPairProvider.hpp"
#include "QuICC/Communicators/Converters/IIndexConv.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of the interface for a data converter.
    */
   class IConverter
   {
      public:
         /// Typedef for forward datatype
         typedef DynamicPairProvider::FwdProvider::RealStorage::DataType RealFwdData;

         /// Typedef for forward datatype
         typedef DynamicPairProvider::FwdProvider::ComplexStorage::DataType ComplexFwdData;

         /// Typdef variant data type for forward data
         typedef std::variant<RealFwdData *, ComplexFwdData *> VariantFwdPointer;

         /// Typedef for backward datatype
         typedef DynamicPairProvider::BwdProvider::RealStorage::DataType RealBwdData;

         /// Typedef for backward datatype
         typedef DynamicPairProvider::BwdProvider::ComplexStorage::DataType ComplexBwdData;

         /// Typdef variant data type for forward data
         typedef std::variant<RealBwdData *, ComplexBwdData *> VariantBwdPointer;

         /**
          * @brief Constructor
          */
         IConverter();

         /**
          * @brief Destructor
          */
         virtual ~IConverter();

         /**
          * @brief Set up the converter
          */
         virtual void setup() = 0;

         /**
          * @brief Convert real data from Fwd to Bwd
          */
         virtual void convertFwd(const VariantFwdPointer &pIn, DynamicPairProvider &storage);

         /**
          * @brief Convert real data from Fwd to Bwd
          */
         virtual void convertFwd(const RealFwdData &in, DynamicPairProvider &storage) = 0;

         /**
          * @brief Convert complex data from Fwd to Bwd
          */
         virtual void convertFwd(const ComplexFwdData &in, DynamicPairProvider &storage) = 0;

         /**
          * @brief Convert real data from Bwd to Fwd
          */
         virtual void convertBwd(const VariantBwdPointer &pIn, DynamicPairProvider &storage);

         /**
          * @brief Convert real data from Bwd to Fwd
          */
         virtual void convertBwd(const RealBwdData &in, DynamicPairProvider &storage) = 0;

         /**
          * @brief Convert complex data from Bwd to Fwd
          */
         virtual void convertBwd(const ComplexBwdData &in, DynamicPairProvider &storage) = 0;

         /**
          * @brief Get the converted real data from Bwd to Fwd
          */
         virtual void getFwd(VariantFwdPointer& pOut, DynamicPairProvider &storage);

         /**
          * @brief Get the converted real data from Bwd to Fwd
          */
         virtual void getFwd(RealFwdData *& pOut, DynamicPairProvider &storage) = 0;

         /**
          * @brief Get the converted complex data from Bwd to Fwd
          */
         virtual void getFwd(ComplexFwdData *& pOut, DynamicPairProvider &storage) = 0;

         /**
          * @brief Get the converted real data from Fwd to Bwd
          */
         virtual void getBwd(VariantBwdPointer & pOut, DynamicPairProvider &storage);

         /**
          * @brief Get the converted real data from Fwd to Bwd
          */
         virtual void getBwd(RealBwdData *& pOut, DynamicPairProvider &storage) = 0;

         /**
          * @brief Get the converted complex data from Fwd to Bwd
          */
         virtual void getBwd(ComplexBwdData *& pOut, DynamicPairProvider &storage) = 0;

         /**
          * @brief Setup upcoming communication
          *
          * @param packs Number of packets in communication packing
          */
         virtual void setupCommunication(const int packs, const TransformDirection::Id direction) = 0;

         /**
          * @brief Start persistent send for forward transform
          */
         virtual void initiateForwardSend() = 0;

         /**
          * @brief Post persistent receive for forward transform
          */
         virtual void prepareForwardReceive() = 0;

         /**
          * @brief Start persistent send for backward transform
          */
         virtual void initiateBackwardSend() = 0;

         /**
          * @brief Post persistent receive for backward transform
          */
         virtual void prepareBackwardReceive() = 0;

         /**
          * @brief Set index converter
          */
         void setIndexConverter(std::shared_ptr<IIndexConv> spConv);

         /**
          * @brief Get index converter
          */
         const IIndexConv& idxConv();

         /**
         * @brief Do storage profiling
         */
         virtual void profileStorage() const = 0;

      protected:
         /**
          * @brief Dimensions of data
          */
         int mDimensions;

      private:
         /**
          * @brief Index converter
          */
         std::shared_ptr<IIndexConv>  mspIdxConv;
   };

   inline const IIndexConv& IConverter::idxConv()
   {
      return *this->mspIdxConv;
   }

   inline void IConverter::convertFwd(const VariantFwdPointer &pIn, DynamicPairProvider &storage)
   {
      std::visit([&](auto&& arg){this->convertFwd(*arg, storage);}, pIn);
   }

   inline void IConverter::convertBwd(const VariantBwdPointer &pIn, DynamicPairProvider &storage)
   {
      std::visit([&](auto&& arg){this->convertBwd(*arg, storage);}, pIn);
   }

   inline void IConverter::getFwd(VariantFwdPointer& pOut, DynamicPairProvider &storage)
   {
      std::visit([&](auto&& arg){this->getFwd(arg, storage);}, pOut);
   }

   inline void IConverter::getBwd(VariantBwdPointer & pOut, DynamicPairProvider &storage)
   {
      std::visit([&](auto&& arg){this->getBwd(arg, storage);}, pOut);
   }

}
}

#endif // QUICC_PARALLEL_ICONVERTER_HPP
