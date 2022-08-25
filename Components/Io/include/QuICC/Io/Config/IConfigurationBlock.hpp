/**
 * @file IConfigurationBlock.hpp
 * @brief Implementation of the base for a configuration file
 */

#ifndef QUICC_IO_CONFIG_ICONFIGURATIONBLOCK_HPP
#define QUICC_IO_CONFIG_ICONFIGURATIONBLOCK_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <map>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Io/Config/IConfigurationNode.hpp"

namespace QuICC {

namespace Io {

namespace Config {

   /**
    * @brief Implementation of the base for a configuration file
    */
   class IConfigurationBlock
   {
      public:
         /// Typedef for node storage
         typedef std::map<std::size_t, SharedIConfigurationNode> NodeMap;

         /**
          * @brief Constructor
          *
          */
         IConfigurationBlock(const std::string& tag);

         /**
          * @brief Destructor
          */
         virtual ~IConfigurationBlock();

         /**
          * @brief Output run information
          */
         void printInfo() const;

         /**
          * @brief Get name of the framework XML tag
          */
         const std::string& tag() const;

         /**
          * @brief Number of nodes
          */
         size_t size() const;

         /**
          * @brief Get iterator range
          */
         std::pair<NodeMap::iterator,NodeMap::iterator> range();

         /**
          * @brief Get const iterator range
          */
         std::pair<NodeMap::const_iterator,NodeMap::const_iterator> crange() const;

         /**
          * @brief Get truncation part
          */
         SharedCIConfigurationNode spNode(const NodeMap::key_type id) const;

         /**
          * @brief Get truncation part
          */
         const IConfigurationNode& node(const NodeMap::key_type id) const;

         /**
          * @brief Get physical part
          */
         SharedIConfigurationNode rspNode(const NodeMap::key_type id);

         /**
          * @brief Gather parameters into arrays
          */
         void gatherParameters(std::vector<int>& iData, std::vector<MHDFloat>& fData, std::vector<std::string>& sData);

         /**
          * @brief Scatter parameters
          */
         void scatterParameters(int& iIdx, int& fIdx, int& sIdx, const std::vector<int>& iData, const std::vector<MHDFloat>& fData, const std::vector<std::string>& sData);

      protected:
         /**
          * @brief Add configuration part to simulation configuration
          *
          * @param id      ID of configuration node
          * @param spNode  Configuration node to add
          */
         void addNode(const NodeMap::key_type id, SharedIConfigurationNode spNode);

         /**
          * @brief Framework part of configuration file
          */
         NodeMap  mBlock;

      private:
         /**
          * @brief XML parent tag name
          */
         std::string mTag;
   };

   /// Typedef for a shared pointer of a configuration file block
   typedef std::shared_ptr<IConfigurationBlock> SharedIConfigurationBlock;

   /// Typedef for a shared pointer of a configuration file node
   typedef std::shared_ptr<const IConfigurationBlock> SharedCIConfigurationBlock;

}
}
}

#endif // QUICC_IO_CONFIG_ICONFIGURATIONBLOCK_HPP
