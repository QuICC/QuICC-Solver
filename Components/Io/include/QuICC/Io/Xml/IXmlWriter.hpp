/**
 * @file IXmlWriter.hpp
 * @brief Implementation of an XML writer
 */

#ifndef QUICC_IO_XML_IXMLWRITER_HPP
#define QUICC_IO_XML_IXMLWRITER_HPP

// Configuration includes
//

// System includes
//
#include <fstream>
#include <sstream>
#include <memory>

// External includes
//
#include "rapidxml.hpp"

// Project includes
//
#include "QuICC/Io/Xml/XmlFile.hpp"

namespace QuICC {

namespace Io {

namespace Xml {

   /**
    * @brief Simpler reader class for XML data (based on rapidXML)
    */
   class IXmlWriter: public XmlFile
   {
      public:
         /**
         * @brief Constructor
         *
         * @param name File name
         * @param ext Extension of the file
         * @param header Header of the file
         * @param type Type of the file
         * @param version Version string of the file
         * @param root XML Root of the file
         */
         IXmlWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, std::string root);

         /**
         * @brief Destructor
         */
         virtual ~IXmlWriter();

         /**
          * @brief Initialise the file
          */
         virtual void init();

         /**
          * @brief Read the content
          */

         virtual void write();

         /**
          * @brief Finalise the file
          */
         virtual void finalize();

      protected:
         /**
          * @brief Handle to the file
          */
         std::ofstream mFile;

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Close the file
          */
         void close();

         /**
          * @brief templated function to write value to XML tree
          *
          * @param val  Value to write
          * @param base XML tree node
          * @param tag  XML tag to write value to
          *
          * \tparam T Type of the value to write
          */
         template <typename T> void writeValue(const T& val, rapidxml::xml_node<>* base, const std::string& tag);

         /**
          * @brief Operation to perform just before writing data
          */
         void preWrite();

         /**
          * @brief Operation to perform just after writing data
          */
         void postWrite();

      private:
   };

   template <typename T>  void IXmlWriter::writeValue(const T& val, rapidxml::xml_node<>* base, const std::string& tag)
   {
      // Pointer to node
      rapidxml::xml_node<>* node;

      // Create node
      node = this->mXML.allocate_node(rapidxml::node_element, tag.c_str());
      base->append_node(node);

      // Create string stream to do type conversion
      std::ostringstream   oss;
      std::string tmp;

      // Convert value into string
      oss << val;

      // Set node value
      char* cstr = this->mXML.allocate_string(oss.str().c_str());
      node->value(cstr);
   }

   /// Typedef for a shared pointer of a IXmlWriter
   typedef std::shared_ptr<IXmlWriter> SharedIXmlWriter;

}
}
}

#endif // QUICC_IO_XML_IXMLWRITER_HPP
