/**
 * @file SerialConverter.cpp
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
#include "QuICC/Communicators/Converters/SerialConverter.hpp"

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Parallel {

   SerialConverter::SerialConverter()
   {
   }

   SerialConverter::~SerialConverter()
   {
   }

   void SerialConverter::getFwd(RealFwdData *& pData, DynamicPairProvider  &storage)
   {
      this->getFwdImpl(pData, storage);
   }

   void SerialConverter::getFwd(ComplexFwdData *& pData, DynamicPairProvider  &storage)
   {
      this->getFwdImpl(pData, storage);
   }

   void SerialConverter::getBwd(RealBwdData *& pData, DynamicPairProvider  &storage)
   {
      this->getBwdImpl(pData, storage);
   }

   void SerialConverter::getBwd(ComplexBwdData *& pData, DynamicPairProvider  &storage)
   {
      this->getBwdImpl(pData, storage);
   }

   void SerialConverter::setup()
   {
   }

   void SerialConverter::setupCommunication(const int packs, const TransformDirection::Id direction)
   {
   }

   void SerialConverter::initiateForwardSend()
   {
   }

   void SerialConverter::prepareForwardReceive()
   {
   }

   void SerialConverter::initiateBackwardSend()
   {
   }

   void SerialConverter::prepareBackwardReceive()
   {
   }

   void SerialConverter::init(SharedResolution spRes, const Dimensions::Transform::Id id, std::shared_ptr<IIndexConv> spIdxConv)
   {
      // Set dimensions
      this->mDimensions = spRes->sim().ss().dimension();

      // Store the shared pointer to the transform resolution
      this->mspTRes = spRes->cpu()->dim(id);

      // Set index converter
      this->setIndexConverter(spIdxConv);
   }

   void SerialConverter::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
#endif // QUICC_STORAGEPROFILE
   }

}
}
