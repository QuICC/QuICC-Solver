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
#include "QuICC/Communicators/Converters/PassthroughIndexConv.hpp"

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
      Parallel::PassthroughIndexConv inConv;
      this->processBwdImpl<RealFwdData,Dimensions::Data::DATF1D>(in, storage, this->idxConv(), inConv);
   }

   void SerialConverterBase::convertFwd(const ComplexFwdData& in, DynamicPairProvider& storage)
   {
      Parallel::PassthroughIndexConv inConv;
      this->processBwdImpl<ComplexFwdData,Dimensions::Data::DATF1D>(in, storage, this->idxConv(), inConv);
   }

   void SerialConverterBase::convertBwd(const RealBwdData& in, DynamicPairProvider& storage)
   {
      Parallel::PassthroughIndexConv outConv;
      this->processFwdImpl<RealBwdData,Dimensions::Data::DATF1D>(in, storage, outConv, this->idxConv());
   }

   void SerialConverterBase::convertBwd(const ComplexBwdData& in, DynamicPairProvider& storage)
   {
      Parallel::PassthroughIndexConv outConv;
      this->processFwdImpl<ComplexBwdData,Dimensions::Data::DATF1D>(in, storage, outConv, this->idxConv());
   }

}
}
