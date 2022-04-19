/**
 * @file IConverter.cpp
 * @brief Source of the interface for a converter
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Communicators/Converters/IConverter.hpp"

// Project includes
//

namespace QuICC {

namespace Parallel {

   IConverter::IConverter()
      : mDimensions(-1)
   {
   }

   IConverter::~IConverter()
   {
   }

   void IConverter::setIndexConverter(std::shared_ptr<IIndexConv> spConv)
   {
      this->mspIdxConv = spConv;
   }

}
}
