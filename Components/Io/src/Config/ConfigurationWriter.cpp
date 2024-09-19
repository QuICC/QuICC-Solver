/**
 * @file ConfigurationWriter.cpp
 * @brief Source of the implementation of the XML parameters file reader
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
#include "QuICC/Io/Config/ConfigurationWriter.hpp"

// Project includes
//
#include "Environment/QuICCEnv.hpp"

namespace QuICC {

namespace Io {

namespace Config {

   const std::string ConfigurationWriter::NAME = "parameters_template_";

   ConfigurationWriter::ConfigurationWriter(const int dim, const std::vector<bool>& isPeriodicBox, const std::string& type)
      : IConfigurationFile<Io::Xml::IXmlWriter>(dim, isPeriodicBox, ConfigurationWriter::NAME+type, type)
   {
   }

   ConfigurationWriter::~ConfigurationWriter()
   {
   }

   void ConfigurationWriter::write()
   {
      // Do pre writing processing
      this->preWrite();

      // Check if the framework allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Write framework block
         this->writeBlock(this->mspFramework);

         // Write setup block
         this->writeBlock(this->mspSetup);

         // Write setup block
         this->writeBlock(this->mspSimulation);

         // Write setup block
         this->writeBlock(this->mspModel);
      }

      // Write xml content to file
      this->mFile << this->mXML;

      // Share Parameters with other CPUs (if applicable)
      this->spreadParameters();

      // Do post writing processing
      this->postWrite();
   }

   void ConfigurationWriter::writeBlock(SharedIConfigurationBlock spBlock)
   {
      //
      // Write simulation configuration
      //
      if(spBlock->size() > 0)
      {
         // Create framework XML node
         rapidxml::xml_node<> *pMaster = this->mXML.allocate_node(rapidxml::node_element, spBlock->tag().c_str());
         this->mpRoot->append_node(pMaster);

         // Check if master node exists
         if(pMaster)
         {
            // Define component node pointer
            rapidxml::xml_node<> *pComponent;

            // Iterate over all master components
            auto itRange = spBlock->crange();
            for(auto itM = itRange.first; itM != itRange.second; itM++)
            {
               if(itM->second->size() > 0)
               {
                  // Create component node
                  pComponent = this->mXML.allocate_node(rapidxml::node_element, itM->second->parent().c_str());
                  pMaster->append_node(pComponent);

                  // Check if component node exists
                  if(pComponent)
                  {
                     // Create integer component iterator
                     auto iRange = itM->second->iTags().crange();
                     // Iterate over all component entries
                     for(auto itIC = iRange.first; itIC != iRange.second; itIC++)
                     {
                        // Create entry value
                        this->writeValue(itIC->second, pComponent, itIC->first);
                     }

                     // Create float component iterator
                     auto fRange = itM->second->fTags().crange();
                     // Iterate over all component entries
                     for(auto itIF = fRange.first; itIF != fRange.second; itIF++)
                     {
                        // Create entry value
                        this->writeValue(itIF->second, pComponent, itIF->first);
                     }

                     // Create float component iterator
                     auto sRange = itM->second->sTags().crange();
                     // Iterate over all component entries
                     for(auto itIS = sRange.first; itIS != sRange.second; itIS++)
                     {
                        // Create entry value
                        this->writeValue(itIS->second, pComponent, itIS->first);
                     }

                     // Check if component data is correct
                     itM->second->checkData();
                  } else
                  {
                     throw std::logic_error("Couldn't find " + itM->second->parent() + " component tag!");
                  }
               }
            }
         } else
         {
            throw std::logic_error("Couldn't find " + spBlock->tag() + " master tag!");
         }
      }
   }
}
}
}
