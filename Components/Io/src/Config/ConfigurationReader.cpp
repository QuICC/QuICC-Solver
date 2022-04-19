/** 
 * @file ConfigurationReader.cpp
 * @brief Source of the implementation of the XML parameters file reader
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
#include "QuICC/Io/Config/ConfigurationReader.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"

namespace QuICC {

namespace Io {

namespace Config {

   const std::string ConfigurationReader::NAME = "parameters";

   ConfigurationReader::ConfigurationReader(const int dim, const std::vector<bool>& isPeriodicBox, const std::string& type)
      : IConfigurationFile<Io::Xml::IXmlReader>(dim, isPeriodicBox, ConfigurationReader::NAME, type)
   {
   }

   ConfigurationReader::~ConfigurationReader()
   {
   }

   void ConfigurationReader::read()
   {
      // Check if the framework allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Read framework configuration
         this->readBlock(this->mspFramework, true);

         // Read setup configuration
         this->readBlock(this->mspSetup, false);

         // Read simulation configuration
         this->readBlock(this->mspSimulation, true);

         // Read model configuration
         this->readBlock(this->mspModel, false);
      }

      // Share Parameters with other CPUs (if applicable)
      this->spreadParameters();
   }

   void ConfigurationReader::readBlock(SharedIConfigurationBlock spBlock, const bool required)
   {
      if(spBlock->size() > 0)
      {
         // Get master pointer to framework XML code
         rapidxml::xml_node<> *pMaster = this->mXML.first_node(spBlock->tag().c_str());

         // Check if master node exists
         if(pMaster)
         {
            // Define component node pointer
            rapidxml::xml_node<> *pComponent;

            // Iterate over all master components
            auto range = spBlock->range();
            for(auto itM = range.first; itM != range.second; itM++)
            {
               if(itM->second->size() > 0)
               {
                  // Get pointer to component node
                  pComponent = pMaster->first_node(itM->second->parent().c_str());

                  // Check if component node exists
                  if(pComponent)
                  {
                     // Create integer component iterator
                     auto iRange = itM->second->iTags().crange();
                     // Iterate over all component entries
                     for(auto itIC = iRange.first; itIC != iRange.second; itIC++)
                     {
                        // Read entry value from file
                        int val;
                        this->readValue(val, pComponent, itIC->first);

                        // Store value
                        itM->second->iTags().setValue(itIC->first, val);
                     }

                     // Create float component iterator
                     auto fRange = itM->second->fTags().crange();
                     // Iterate over all component entries
                     for(auto itIF = fRange.first; itIF != fRange.second; itIF++)
                     {
                        // Read entry value from file
                        MHDFloat val;
                        this->readValue(val, pComponent, itIF->first);

                        // Store value
                        itM->second->fTags().setValue(itIF->first, val);
                     }

                     // Check if component data is correct
                     itM->second->checkData();
                  } else
                  {
                     throw std::logic_error("Couldn't find " + itM->second->parent() + " component tag from configuration file!");
                  }
               }
            }
         } else
         {
            if(required)
            {
               throw std::logic_error("Couldn't find " + spBlock->tag() + " master tag from configuration file!");
            }
         }
      }
   }
}
}
}
