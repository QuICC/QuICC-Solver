/** 
 * @file IConfigurationBlock.cpp
 * @brief Source of the implementation of the base of the configuration file block
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Config/IConfigurationBlock.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Tools/Formatter.hpp"

namespace QuICC {

namespace Io {

namespace Config {

   IConfigurationBlock::IConfigurationBlock(const std::string& tag)
      : mTag(tag)
   {
   }

   IConfigurationBlock::~IConfigurationBlock()
   {
   }

   size_t IConfigurationBlock::size() const
   {
      return this->mBlock.size();
   }

   std::pair<IConfigurationBlock::NodeMap::iterator,IConfigurationBlock::NodeMap::iterator> IConfigurationBlock::range()
   {
      auto range = std::make_pair(this->mBlock.begin(), this->mBlock.end());
      return range;
   }

   std::pair<IConfigurationBlock::NodeMap::const_iterator,IConfigurationBlock::NodeMap::const_iterator> IConfigurationBlock::crange() const
   {
      auto range = std::make_pair(this->mBlock.cbegin(), this->mBlock.cend());
      return range;
   }

   SharedCIConfigurationNode IConfigurationBlock::spNode(const NodeMap::key_type id) const
   {
      // Make sure initialisation was correct
      assert(this->mBlock.find(id) != this->mBlock.end());

      return this->mBlock.find(id)->second;
   }

   const IConfigurationNode& IConfigurationBlock::node(const NodeMap::key_type id) const
   {
      // Make sure initialisation was correct
      assert(this->mBlock.find(id) != this->mBlock.end());

      return *this->mBlock.find(id)->second;
   }

   SharedIConfigurationNode IConfigurationBlock::rspNode(const NodeMap::key_type id)
   {
      // Make sure initialisation was correct
      assert(this->mBlock.find(id) != this->mBlock.end());

      return this->mBlock.find(id)->second;
   }

   const std::string& IConfigurationBlock::tag() const
   {
      return this->mTag;
   }

   void IConfigurationBlock::addNode(const NodeMap::key_type id, SharedIConfigurationNode spNode)
   {
      this->mBlock.insert(std::make_pair(id, spNode));
   }

   void IConfigurationBlock::printInfo() const
   {
      // Check if the framework allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Create output header
         Tools::Formatter::printNewline(std::cout);
         Tools::Formatter::printLine(std::cout, '-');
         Tools::Formatter::printCentered(std::cout, this->mTag, '*');
         Tools::Formatter::printLine(std::cout, '-');
         Tools::Formatter::printNewline(std::cout);

         // Iterate over all nodes
         for(auto it = this->mBlock.cbegin(); it != this->mBlock.cend(); it++)
         {
            it->second->printInfo();
         }
      }
   }

   void IConfigurationBlock::gatherParameters(std::vector<int>& iData, std::vector<MHDFloat>& fData)
   {
      // Iterate over all master components
      for(auto citB = this->mBlock.cbegin(); citB != this->mBlock.cend(); citB++)
      {
         // Iterate over all component entries
         auto iRange = citB->second->iTags().crange();
         for(auto citIC = iRange.first; citIC != iRange.second; citIC++)
         {
            iData.push_back(citIC->second);
         }

         // Iterate over all component entries
         auto fRange = citB->second->fTags().crange();
         for(auto citIF = fRange.first; citIF != fRange.second; citIF++)
         {
            fData.push_back(citIF->second);
         }
      }
   }

   void IConfigurationBlock::scatterParameters(int& iIdx, int& fIdx, const std::vector<int>& iData, const std::vector<MHDFloat>& fData)
   {
      // Iterate over all master components
      for(auto itB = this->mBlock.begin(); itB != this->mBlock.end(); itB++)
      {
         // Create integer component iterator
         auto iRange = itB->second->iTags().crange();
         // Iterate over all component entries
         for(auto itIC = iRange.first; itIC != iRange.second; itIC++)
         {
            itB->second->iTags().setValue(itIC->first, iData.at(iIdx));
            // Increment integer index
            iIdx++;
         }

         // Create float component iterator
         auto fRange = itB->second->fTags().crange();
         // Iterate over all component entries
         for(auto itIF = fRange.first; itIF != fRange.second; itIF++)
         {
            itB->second->fTags().setValue(itIF->first, fData.at(fIdx));
            // Increment float index
            fIdx++;
         }
      }
   }


}
}
}
