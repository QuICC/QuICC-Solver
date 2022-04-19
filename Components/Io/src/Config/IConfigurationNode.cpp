/** 
 * @file IConfigurationNode.cpp
 * @brief Source of the implementation of the base of the configuration file nodes
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Config/IConfigurationNode.hpp"

// Project includes
//
#include "QuICC/Tools/Formatter.hpp"

namespace QuICC {

namespace Io {

namespace Config {

   IConfigurationNode::IConfigurationNode(const std::string& parent)
      : mParent(parent)
   {
   }

   IConfigurationNode::~IConfigurationNode()
   {
   }

   const std::string& IConfigurationNode::parent() const
   {
      return this->mParent;
   }

   size_t IConfigurationNode::size() const
   {
      return this->mIData.size() + this->mFData.size() + mSData.size();
   }

   IConfigurationTag<int>& IConfigurationNode::iTags()
   {
      return this->mIData;
   }

   const IConfigurationTag<int>& IConfigurationNode::iTags() const
   {
      return this->mIData;
   }

   IConfigurationTag<MHDFloat>& IConfigurationNode::fTags()
   {
      return this->mFData;
   }

   const IConfigurationTag<MHDFloat>& IConfigurationNode::fTags() const
   {
      return this->mFData;
   }

   IConfigurationTag<std::string>& IConfigurationNode::sTags()
   {
      return this->mSData;
   }

   const IConfigurationTag<std::string>& IConfigurationNode::sTags() const
   {
      return this->mSData;
   }

   void IConfigurationNode::printInfo() const
   {
      // Create header
      Tools::Formatter::printLine(std::cout, '-');
      Tools::Formatter::printCentered(std::cout, this->parent(), '*');
      Tools::Formatter::printLine(std::cout, '-');

      std::stringstream oss;

      this->mIData.printInfo();

      this->mFData.printInfo();

      this->mSData.printInfo();

      Tools::Formatter::printLine(std::cout, '-');
      Tools::Formatter::printNewline(std::cout);
   }

}
}
}
