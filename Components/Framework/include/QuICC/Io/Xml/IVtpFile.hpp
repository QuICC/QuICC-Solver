/**
 * @file IVtpFile.hpp
 * @brief Implementation of the base for a VTK PolyData XML format file
 */

#ifndef QUICC_IO_XML_IVTPFILE_HPP
#define QUICC_IO_XML_IVTPFILE_HPP

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
    * @brief Implementation of the base for a VTK PolyData XML format file
    */
   template <typename TBase> class IVtpFile: public TBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name Name of the file
          */
         IVtpFile(const std::string& name);

         /**
          * @brief Destructor
          */
         virtual ~IVtpFile();

         /**
          * @brief Initialise the file
          */
         virtual void init();

      protected:
         /**
          * @brief EXTENSION of the VTK XML file
          */
         static const std::string EXTENSION;

         /**
          * @brief HEADER of the VTK XML file
          */
         static const std::string HEADER;

         /**
          * @brief TYPE of the VTK XML file
          */
         static const std::string TYPE;

         /**
          * @brief VERSION of the VTK XML file
          */
         static const std::string VERSION;

      private:
   };

   template <typename TBase> const std::string IVtpFile<TBase>::EXTENSION = ".vtp";

   template <typename TBase> const std::string IVtpFile<TBase>::HEADER = "VTKFile";

   template <typename TBase> const std::string IVtpFile<TBase>::TYPE = "PolyData";

   template <typename TBase> const std::string IVtpFile<TBase>::VERSION = "1.0";

   template <typename TBase> IVtpFile<TBase>::IVtpFile(const std::string& name)
      : TBase(name, IVtpFile<TBase>::EXTENSION, IVtpFile<TBase>::HEADER, IVtpFile<TBase>::TYPE, IVtpFile<TBase>::VERSION)
   {
   }

   template <typename TBase> IVtpFile<TBase>::~IVtpFile()
   {
   }

   template <typename TBase> void IVtpFile<TBase>::init()
   {
      // Node variables
      rapidxml::xml_node<>* pNode;
      rapidxml::xml_node<>* pVtp;

      // Create xml declaration
      pNode = this->mXML.allocate_node(rapidxml::node_declaration);
      pNode->append_attribute(this->mXML.allocate_attribute("version", "1.0"));
      pNode->append_attribute(this->mXML.allocate_attribute("encoding", "utf-8"));
      this->mXML.append_node(pNode);

      // VTKFile node
      pVtp = this->mXML.allocate_node(rapidxml::node_element, this->HEADER.c_str());
      pVtp->append_attribute(this->mXML.allocate_attribute("type", this->TYPE.c_str()));
      pVtp->append_attribute(this->mXML.allocate_attribute("version", this->VERSION.c_str()));
      pVtp->append_attribute(this->mXML.allocate_attribute("byte_order", "LittleEndian"));
      this->mXML.append_node(pVtp);

      // Create PolyData node
      pNode = this->mXML.allocate_node(rapidxml::node_element, this->TYPE.c_str());
      pVtp->append_node(pNode);
   }

}
}
}

#endif // QUICC_IO_XML_IVTPFILE_HPP
