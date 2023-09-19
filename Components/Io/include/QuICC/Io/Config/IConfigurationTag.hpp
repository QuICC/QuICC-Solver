/**
 * @file IConfigurationTag.hpp
 * @brief Implementation of a configuration tag of the configuration file
 */

#ifndef QUICC_IO_CONFIG_ICONFIGURATIONTAG_HPP
#define QUICC_IO_CONFIG_ICONFIGURATIONTAG_HPP

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
#include "QuICC/Tools/Formatter.hpp"

namespace QuICC {

namespace Io {

namespace Config {

   /**
    * @brief Implementation of a configuration tag of the configuration file
    */
   template <typename T> class IConfigurationTag
   {
      public:
         /// Map type
         typedef std::map<std::string,T> MapType;
         /// Map const iterator range
         typedef std::pair<typename MapType::const_iterator,typename MapType::const_iterator> ConstMapRangeType;

         /**
          * @brief Constructor
          */
         explicit IConfigurationTag();

         /**
          * @brief Destructor
          */
         virtual ~IConfigurationTag();

         /**
          * @brief Get size
          */
         size_t size() const;

         /**
          * @brief Get integer value by name
          *
          * @param name Name of the value
          */
         T value(std::string name) const;

         /**
          * @brief Set integer value by name
          *
          * @param name    Name of the value
          * @param value   New value
          */
         void setValue(std::string name, T value);

         /**
          * @brief Get integer name to value map
          */
         const MapType& map() const;

         /**
          * @brief Get integer name to value map
          */
         MapType& rMap();

         /**
          * @brief Get integer name to value map iterator range
          */
         ConstMapRangeType crange() const;

         /**
          * @brief Output run information
          */
         void printInfo() const;

         /**
          * @brief Add XML tag with integer data
          */
         void addTag(const std::string& name, const T value);

      protected:

      private:
         /**
          * @brief XML tags of the integer data
          */
         MapType mData;
   };

   template <typename T> IConfigurationTag<T>::IConfigurationTag()
   {
   }

   template <typename T> IConfigurationTag<T>::~IConfigurationTag()
   {
   }

   template <typename T> size_t IConfigurationTag<T>::size() const
   {
      return this->mData.size();
   }

   template <typename T> T IConfigurationTag<T>::value(std::string name) const
   {
      // Make sure initialisation was correct
      assert(this->mData.find(name) != this->mData.end());

      return this->mData.find(name)->second;
   }

   template <typename T> void IConfigurationTag<T>::setValue(std::string name, const T value)
   {
      // Make sure initialisation was correct
      assert(this->mData.find(name) != this->mData.end());

      this->mData.find(name)->second = value;
   }

   template <typename T> const typename IConfigurationTag<T>::MapType& IConfigurationTag<T>::map() const
   {
      return this->mData;
   }

   template <typename T> typename IConfigurationTag<T>::MapType& IConfigurationTag<T>::rMap()
   {
      return this->mData;
   }

   template <typename T> typename IConfigurationTag<T>::ConstMapRangeType IConfigurationTag<T>::crange() const
   {
      return std::make_pair(this->mData.cbegin(),this->mData.cend());
   }

   template <typename T> void IConfigurationTag<T>::printInfo() const
   {
      std::stringstream oss;
      for(auto it = this->mData.cbegin(); it != this->mData.cend(); it++)
      {
         oss << it->first << ": " << it->second;
         Tools::Formatter::printCentered(std::cout, oss.str(), ' ');
         oss.str("");
      }
   }

   template <typename T> void IConfigurationTag<T>::addTag(const std::string& name, T value)
   {
      this->mData.insert(std::make_pair(name, value));
   }

}
}
}

#endif // QUICC_IO_CONFIG_ICONFIGURATIONTAG_HPP
