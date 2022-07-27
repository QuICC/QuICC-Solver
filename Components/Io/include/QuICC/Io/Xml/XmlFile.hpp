/**
 * @file XmlFile.hpp
 * @brief Implementation of a general XML file
 */

#ifndef QUICC_IO_XML_XMLFILE_HPP
#define QUICC_IO_XML_XMLFILE_HPP

// System includes
//
#include <string>

// External includes
//
#include <rapidxml.hpp>

// Project includes
//

namespace QuICC {

namespace Io {

namespace Xml {

   /**
    * @brief Implementation of a general XML file
    */
   class XmlFile
   {
      public:
         /**
         * @brief Constructor
         *
         * @param name     Filename
         * @param ext      File extension
         * @param header   Header string of file
         * @param type     Type string of file
         * @param version  Version string of file
         * @param root     XML root
         */
         XmlFile(std::string name, std::string ext, std::string header, std::string type, std::string version, std::string root);

         /**
         * @brief Destructor
         */
         virtual ~XmlFile();

         /**
          * @brief Get filename
          */
         std::string  filename() const;

      protected:
         /**
          * @brief XML interface
          */
         rapidxml::xml_document<>   mXML;

         /**
          * @brief XML interface
          */
         rapidxml::xml_node<>*   mpRoot;

         /**
          * @brief Get the name
          */
         const std::string&  name() const;

         /**
          * @brief Get extension
          */
         const std::string&  extension() const;

         /**
          * @brief Get file meta tag
          */
         const std::string& fileTag() const;

         /**
          * @brief Get the header tag
          */
         const std::string&  headerTag() const;

         /**
          * @brief Get the header content
          */
         const std::string&  header() const;

         /**
          * @brief Get the type tag
          */
         const std::string&  typeTag() const;

         /**
          * @brief Get the type content
          */
         const std::string&  type() const;

         /**
          * @brief Get the version tag
          */
         const std::string&  versionTag() const;

         /**
          * @brief Get the version content
          */
         const std::string&  version() const;

         /**
          * @brief Get the root content
          */
         const std::string&  root() const;

         /**
          * @brief Reset name
          *
          * @param name New name
          */
         void resetName(std::string name);

      private:
         /**
          * @brief XML tag for meta file information
          */
         static const std::string FILE_TAG;

         /**
          * @brief XML tag for file header information
          */
         static const std::string HEADER_TAG;

         /**
          * @brief XML tag for file type information
          */
         static const std::string TYPE_TAG;

         /**
          * @brief XML tag for file version information
          */
         static const std::string VERSION_TAG;

         /**
          * @brief Name of the file without extension
          */
         std::string mName;

         /**
          * @brief File extension
          */
         std::string mExt;

         /**
          * @brief Header of the file to check compatibility
          */
         std::string mHeader;

         /**
          * @brief Type of the file to check compatibility
          */
         std::string mType;

         /**
          * @brief Version of the file to check compatibility
          */
         std::string mVersion;

         /**
          * @brief Version of the file to check compatibility
          */
         std::string mRoot;
   };

}
}
}

#endif // QUICC_IO_XML_XMLFILE_HPP
