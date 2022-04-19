/** 
 * @file IXmlReader.cpp 
 * @brief Source of the interface to the XML reader
 */

// Configuration includes
//

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Io/Xml/IXmlReader.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"

namespace QuICC {

namespace Io {

namespace Xml {

   IXmlReader::IXmlReader(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : XmlFile(name, ext, header, type, version), mContent()
   {
   }

   IXmlReader::~IXmlReader()
   {
   }

   void IXmlReader::init()
   {
      // Check if the framework allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Open the file
         this->open();

         // Parse the XML tree
         this->parse();

         // Check tree compatibility
         this->checkCompatibility();
      }
   }

   void IXmlReader::open()
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

   void IXmlReader::parse()
   {
      // Check if the framework allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         std::string tmp;
         std::string tmpAll;

         // Read file, line by line until EOF 
         while(this->mFile.good())
         {
            getline(this->mFile, tmp);
            tmpAll += tmp;
         }

         // Create vector of char out of file content. Adding end character
         this->mContent = std::vector<char>(tmpAll.begin(), tmpAll.end());
         this->mContent.push_back('\0');

         // Create XML tree through parser out of char content
         this->mXML.parse<0>(&this->mContent[0]);
      }
   }

   void IXmlReader::checkCompatibility()
   {
      // Check if the framework allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Get pointer to FILEMETA node
         rapidxml::xml_node<> *node = this->mXML.first_node(this->fileTag().c_str());

         // Get pointer to HEADER subnode
         rapidxml::xml_node<> *vnode = node->first_node(this->headerTag().c_str());

         // Check for successful extraction of HEADER node
         if(vnode)
         {
            // Get HEADER node content
            std::string fileHead = vnode->value();

            // Check HEADER compatibility
            if(fileHead == this->header())
            {
               // Get pointer to TYPE subnode
               vnode = node->first_node(this->typeTag().c_str());

               // Check for successful extraction of TYPE node
               if(vnode)
               {
                  // Get TYPE node content
                  std::string fileType = vnode->value();

                  // Check TYPE compatibility
                  if(fileType == this->type())
                  {
                     // Get pointer to VERSION subnode
                     vnode = node->first_node(this->versionTag().c_str());

                     // Check for successful extraction of VERSION node
                     if(vnode)
                     {
                        // Get VERSION node content
                        std::string fileVers = vnode->value();

                        // Check VERSION compatibility
                        if(fileVers == this->version())
                        {
                           //
                           // File is compatible with reader.
                           // Nothing to do
                           //
                        } else
                        {
                           throw std::logic_error("Wrong XML file version!");
                        }
                     } else
                     {
                        throw std::logic_error("Missing XML file version!");
                     }
                  } else
                  {
                     throw std::logic_error("Wrong XML file type!");
                  }
               } else
               {
                  throw std::logic_error("Missing XML file type!");
               }
            } else
            {
               throw std::logic_error("Wrong XML file header!");
            }
         } else
         {
            throw std::logic_error("Missing XML file header!");
         }
      }
   }

   void IXmlReader::finalize()
   {
      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Close file
         this->close();
      }
   }

   void IXmlReader::close()
   {
      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile.close();
      }
   }

}
}
}
