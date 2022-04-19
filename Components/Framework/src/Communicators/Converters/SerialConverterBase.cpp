/**
 * @file SerialConverterBase.cpp
 * @brief Source of the serial data converter base
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Communicators/Converters/SerialConverterBase.hpp"

// Project includes
//

namespace QuICC {

namespace Parallel {

   SerialConverterBase::SerialConverterBase()
      : IConverter()
   {
   }

   SerialConverterBase::~SerialConverterBase()
   {
   }

   void SerialConverterBase::convertFwd(const RealFwdData& in, DynamicPairProvider& storage)
   {
      this->convertFwdImpl(in, storage);
   }

   void SerialConverterBase::convertFwd(const ComplexFwdData& in, DynamicPairProvider& storage)
   {
      this->convertFwdImpl(in, storage);
   }

   void SerialConverterBase::convertBwd(const RealBwdData& in, DynamicPairProvider& storage)
   {
      this->convertBwdImpl(in, storage);
   }

   void SerialConverterBase::convertBwd(const ComplexBwdData& in, DynamicPairProvider& storage)
   {
      this->convertBwdImpl(in, storage);
   }

}
}
