/** 
 * @file XmlFile.cpp
 * @brief Source of the general XML file implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Xml/XmlFile.hpp"

// Project includes
//

namespace QuICC {

namespace Io {

namespace Xml {

   XmlFile::XmlFile(std::string name, std::string ext, std::string header, std::string type, std::string version, std::string root)
      : mXML(), mName(name), mExt(ext), mHeader(header), mType(type), mVersion(version), mRoot(root)
   {
   }

   XmlFile::~XmlFile()
   {
   }

   std::string XmlFile::filename() const
   {
      return this->mName + this->mExt;
   }

   void XmlFile::resetName(std::string name)
   {
      this->mName = name;
   }

   const std::string& XmlFile::name() const
   {
      return this->mName;
   }

   const std::string& XmlFile::extension() const
   {
      return this->mExt;
   }

   const std::string& XmlFile::fileTag() const
   {
      return XmlFile::FILE_TAG;
   }

   const std::string& XmlFile::headerTag() const
   {
      return XmlFile::HEADER_TAG;
   }

   const std::string& XmlFile::header() const
   {
      return this->mHeader;
   }

   const std::string& XmlFile::typeTag() const
   {
      return XmlFile::TYPE_TAG;
   }

   const std::string& XmlFile::type() const
   {
      return this->mType;
   }

   const std::string& XmlFile::versionTag() const
   {
      return XmlFile::VERSION_TAG;
   }

   const std::string& XmlFile::version() const
   {
      return this->mVersion;
   }

   const std::string& XmlFile::root() const
   {
      return this->mRoot;
   }

   const std::string XmlFile::FILE_TAG = "file";

   const std::string XmlFile::HEADER_TAG = "header";

   const std::string XmlFile::TYPE_TAG = "type";

   const std::string XmlFile::VERSION_TAG = "version";
}
}
}
