/**
 * @file GxlWriter.hpp
 * @brief Implementation of the GXL format file writer
 */

#ifndef QUICC_IO_XML_GXLWRITER_HPP
#define QUICC_IO_XML_GXLWRITER_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Typedefs.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Io/Xml/IXmlWriter.hpp"
#include "QuICC/Io/Xml/IGxlFile.hpp"
#include "QuICC/Enums/TransformDirection.hpp"
#include "QuICC/TransformConfigurators/TransformPath.hpp"
#include "QuICC/TransformConfigurators/TransformTree.hpp"

namespace QuICC {

namespace Io {

namespace Xml {

   /**
    * @brief Implementation of the GXL format file writer
    */
   class GxlWriter: public IGxlFile<Io::Xml::IXmlWriter>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name Name of the file
          */
         GxlWriter(const std::string& name);

         /**
          * @brief Destructor
          */
         virtual ~GxlWriter();

         /**
          * @brief Read content of configuration file
          */
         virtual void write();

         /**
          * @brief Create communication graph
          */
         void graphCommunication(const std::map<Dimensions::Transform::Id,std::multimap<int,int> >& structure);

         /**
          * @brief Create transform path graph
          */
         void graphTransformPath(const std::map<std::size_t,std::vector<Transform::TransformPath> >& paths, const TransformDirection::Id dir);

         /**
          * @brief Create transform tree graph
          */
         void graphTransformTree(const std::vector<Transform::TransformTree>& trees, const TransformDirection::Id dir);

      protected:
         /**
          * @brief Create an attr tag in the xml tree
          *
          * @param parent  Parent node to attach to
          * @param name    Name of the attr tag
          * @param value   Value to store in string child
          */
         void createAttr(rapidxml::xml_node<>* parent, const std::string& name, const std::string& value);

         /**
          * @brief Build tree recursively
          */
         void graphTransformTreeEdge(const Transform::TransformTreeEdge& edge, const std::string& root, std::vector<std::string>::const_iterator colorIt, rapidxml::xml_node<> * pGraph, const TransformDirection::Id dir);

      private:
   };

   /// Typedef for a smart pointer of a GxlWriter
   typedef std::shared_ptr<GxlWriter> SharedGxlWriter;

}
}
}

#endif // QUICC_IO_XML_GXLWRITER_HPP
