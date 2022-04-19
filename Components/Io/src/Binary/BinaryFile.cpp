/** 
 * @file BinaryFile.cpp
 * @brief Source of the general binary file implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Binary/BinaryFile.hpp"

// Project includes
//

namespace QuICC {

namespace Io {

namespace Binary {

   BinaryFile::BinaryFile(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : mName(name), mExt(ext), mHeader(header), mType(type), mVersion(version)
   {
   }

   BinaryFile::~BinaryFile()
   {
   }

   std::string BinaryFile::filename() const
   {
      return this->mName + this->mExt;
   }

   void BinaryFile::resetName(std::string name)
   {
      this->mName = name;
   }

   std::string BinaryFile::name() const
   {
      return this->mName;
   }

   std::string BinaryFile::extension() const
   {
      return this->mExt;
   }

   std::string BinaryFile::header() const
   {
      return BinaryFile::FILE_HEADER + this->mHeader;
   }

   std::string BinaryFile::type() const
   {
      return BinaryFile::FILE_TYPE + this->mType;
   }

   std::string BinaryFile::version() const
   {
      return BinaryFile::FILE_VERSION + this->mVersion;
   }

   const std::string BinaryFile::FILE_HEADER = "#";

   const std::string BinaryFile::FILE_TYPE = "#Type: ";

   const std::string BinaryFile::FILE_VERSION = "#Version: ";
}
}
}
