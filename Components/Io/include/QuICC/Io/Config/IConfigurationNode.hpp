/**
 * @file IConfigurationNode.hpp
 * @brief Implementation of a configuration node of the configuration file
 */

#ifndef QUICC_IO_CONFIG_ICONFIGURATIONNODE_HPP
#define QUICC_IO_CONFIG_ICONFIGURATIONNODE_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <map>
#include <memory>

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Io/Config/IConfigurationTag.hpp"

namespace QuICC {

namespace Io {

namespace Config {

   /**
    * @brief Implementation of a configuration node of the configuration file
    */
   class IConfigurationNode
   {
      public:
         /**
          * @brief Constructor
          */
         explicit IConfigurationNode(const std::string& parent);

         /**
          * @brief Destructor
          */
         virtual ~IConfigurationNode();

         /**
          * @brief Get parent node XML tag
          */
         const std::string& parent() const;

         /**
          * @brief Get size
          */
         size_t size() const;

         /**
          * @brief Set integer tags
          */
         IConfigurationTag<int>& iTags();

         /**
          * @brief Get integer tags
          */
         const IConfigurationTag<int>& iTags() const;

         /**
          * @brief Set float tags
          */
         IConfigurationTag<MHDFloat>& fTags();

         /**
          * @brief Get float tags
          */
         const IConfigurationTag<MHDFloat>& fTags() const;

         /**
          * @brief Set string tags
          */
         IConfigurationTag<std::string>& sTags();

         /**
          * @brief Get string tags
          */
         const IConfigurationTag<std::string>& sTags() const;

         /**
          * @brief Check compatibility of data
          */
         virtual void checkData() = 0;

         /**
          * @brief Output run information
          */
         virtual void printInfo() const;

      protected:

      private:
         /**
          * @brief XML parent tag name
          */
         std::string mParent;

         /**
          * @brief XML tags of the integer data
          */
         IConfigurationTag<int> mIData;

         /**
          * @brief XML tags of the float data
          */
         IConfigurationTag<MHDFloat> mFData;

         /**
          * @brief XML tags of the string data
          */
         IConfigurationTag<std::string> mSData;
   };

   /// Typedef for a shared pointer of a configuration file node
   typedef std::shared_ptr<IConfigurationNode> SharedIConfigurationNode;

   /// Typedef for a shared pointer of a configuration file node
   typedef std::shared_ptr<const IConfigurationNode> SharedCIConfigurationNode;

}
}
}

#endif // QUICC_IO_CONFIG_ICONFIGURATIONNODE_HPP
