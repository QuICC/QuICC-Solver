/** 
 * @file AsciiFile.cpp
 * @brief Source of the general ASCII file implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Ascii/AsciiFile.hpp"

// Project includes
//

namespace QuICC {

namespace Io {

namespace Ascii {

   AsciiFile::AsciiFile(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : mName(name), mExt(ext), mHeader(header), mType(type), mVersion(version)
   {
   }

   AsciiFile::~AsciiFile()
   {
   }

   std::string AsciiFile::filename() const
   {
      return this->mName + this->mExt;
   }

   void AsciiFile::resetName(std::string name)
   {
      this->mName = name;
   }

   std::string AsciiFile::name() const
   {
      return this->mName;
   }

   std::string AsciiFile::extension() const
   {
      return this->mExt;
   }

   std::string AsciiFile::header() const
   {
      return AsciiFile::FILE_HEADER + this->mHeader;
   }

   std::string AsciiFile::type() const
   {
      return AsciiFile::FILE_TYPE + this->mType;
   }

   std::string AsciiFile::version() const
   {
      return AsciiFile::FILE_VERSION + this->mVersion;
   }

   const std::string AsciiFile::FILE_HEADER = "#";

   const std::string AsciiFile::FILE_TYPE = "#Type: ";

   const std::string AsciiFile::FILE_VERSION = "#Version: ";
}
}
}
