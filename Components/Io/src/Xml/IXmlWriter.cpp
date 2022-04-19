/** 
 * @file IXmlWriter.cpp 
 * @brief Source of the implementation of the XML writer
 */

// Configuration includes
//

// System includes
//
#include <stdexcept>

// External includes
//
#include <rapidxml_print.hpp>

// Class include
//
#include "QuICC/Io/Xml/IXmlWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"

namespace QuICC {

namespace Io {

namespace Xml {

   IXmlWriter::IXmlWriter(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : XmlFile(name, ext, header, type, version)
   {
   }

   IXmlWriter::~IXmlWriter()
   {
   }

   void IXmlWriter::init()
   {
      // Node variables
      rapidxml::xml_node<>* node;
      rapidxml::xml_node<>* child;

      // Create xml declaration
      node = this->mXML.allocate_node(rapidxml::node_declaration);
      std::string sver = "version";
      std::string svern = "1.0";
      node->append_attribute(this->mXML.allocate_attribute(this->mXML.allocate_string(sver.c_str()), this->mXML.allocate_string(svern.c_str())));
      std::string senc = "encoding";
      std::string senct = "utf-8";
      node->append_attribute(this->mXML.allocate_attribute(this->mXML.allocate_string(senc.c_str()), this->mXML.allocate_string(senct.c_str())));
      this->mXML.append_node(node);
      
      // FILEMETA node
      node = this->mXML.allocate_node(rapidxml::node_element, this->fileTag().c_str());
      this->mXML.append_node(node);
      
      // HEADER node
      child = this->mXML.allocate_node(rapidxml::node_element, this->headerTag().c_str());
      child->value(this->header().c_str());
      node->append_node(child);
      
      // TYPE node
      child = this->mXML.allocate_node(rapidxml::node_element, this->typeTag().c_str());
      child->value(this->type().c_str());
      node->append_node(child);
      
      // VERSION node
      child = this->mXML.allocate_node(rapidxml::node_element, this->versionTag().c_str());
      child->value(this->version().c_str());
      node->append_node(child);
   }

   void IXmlWriter::open()
   {
      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Get handle to file
         this->mFile.open(this->filename().c_str());

         // Check that opening was a success
         if(! this->mFile.is_open())
         {
            throw std::logic_error("Couldn't open XML file " + this->filename() + "!");
         }
      }
   }

   void IXmlWriter::finalize()
   {
   }

   void IXmlWriter::close()
   {
      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile.close();
      }
   }

   void IXmlWriter::write()
   {
      // Do pre writing processing
      this->preWrite();

      // Check if the framework allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Write xml content to file
         this->mFile << this->mXML;
      }

      // Do post writing processing
      this->postWrite();
   }

   void IXmlWriter::preWrite()
   {
      // Check if the framework allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Create the file
         this->open();
      }
   }

   void IXmlWriter::postWrite()
   {
      // Check if the framework allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Close the file
         this->close();
      }
   }

}
}
}
