/**
 * @file IGxlFile.hpp
 * @brief Implementation of the base for a GXL format file
 */

#ifndef QUICC_IO_XML_IGXLFILE_HPP
#define QUICC_IO_XML_IGXLFILE_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//
#include <rapidxml_print.hpp>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Tools/Formatter.hpp"

namespace QuICC {

namespace Io {

namespace Xml {

   /**
    * @brief Implementation of the base for a GXL format file
    */
   template <typename TBase> class IGxlFile: public TBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name Name of the file
          */
         IGxlFile(const std::string& name);

         /**
          * @brief Destructor
          */
         virtual ~IGxlFile();

         /**
          * @brief Initialise the file
          */
         virtual void init();

      protected:
         /**
          * @brief GXL tag
          */
         static const std::string GXLTAG;

      private:
         /**
          * @brief EXTENSION of the GXL file
          */
         static const std::string EXTENSION;

         /**
          * @brief HEADER of the GXL file
          */
         static const std::string HEADER;

         /**
          * @brief TYPE of the GXL file
          */
         static const std::string TYPE;

         /**
          * @brief VERSION of the GXL file
          */
         static const std::string VERSION;
   };

   template <typename TBase> const std::string IGxlFile<TBase>::EXTENSION = ".gxl";

   template <typename TBase> const std::string IGxlFile<TBase>::HEADER = "GXLFile";

   template <typename TBase> const std::string IGxlFile<TBase>::TYPE = "GraphFile";

   template <typename TBase> const std::string IGxlFile<TBase>::VERSION = "1.0";

   template <typename TBase> const std::string IGxlFile<TBase>::GXLTAG = "gxl";

   template <typename TBase> IGxlFile<TBase>::IGxlFile(const std::string& name)
      : TBase(name, IGxlFile<TBase>::EXTENSION, IGxlFile<TBase>::HEADER, IGxlFile<TBase>::TYPE, IGxlFile<TBase>::VERSION, IGxlFile<TBase>::GXLTAG)
   {
   }

   template <typename TBase> IGxlFile<TBase>::~IGxlFile()
   {
   }

   template <typename TBase> void IGxlFile<TBase>::init()
   {
      // Node variables
      rapidxml::xml_node<>* node;

      // Create xml declaration
      node = this->mXML.allocate_node(rapidxml::node_declaration);
      node->append_attribute(this->mXML.allocate_attribute("version", "1.0"));
      node->append_attribute(this->mXML.allocate_attribute("encoding", "utf-8"));
      this->mXML.append_node(node);

      // DOCTYPE node
      node = this->mXML.allocate_node(rapidxml::node_doctype);
      node->value("gxl SYSTEM \"http://www.gupro.de/GXL/gxl-1.0.dtd\"");
      this->mXML.append_node(node);

      // gxl node
      this->mpRoot = this->mXML.allocate_node(rapidxml::node_element, this->GXLTAG.c_str());
      this->mpRoot->append_attribute(this->mXML.allocate_attribute("xmlns:xlink", "http://www.w3.org/1999/xlink"));
      this->mXML.append_node(this->mpRoot);
   }

}
}
}

#endif // QUICC_IO_XML_IGXLFILE_HPP
