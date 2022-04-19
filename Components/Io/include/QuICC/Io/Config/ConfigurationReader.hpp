/** 
 * @file ConfigurationReader.hpp 
 * @brief Implementation of the XML configuration file reader
 */

#ifndef QUICC_IO_CONFIG_CONFIGURATIONREADER_HPP
#define QUICC_IO_CONFIG_CONFIGURATIONREADER_HPP

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
#include "QuICC/Typedefs.hpp"
#include "QuICC/Io/Xml/IXmlReader.hpp"
#include "QuICC/Io/Config/IConfigurationFile.hpp"

namespace QuICC {

namespace Io {

namespace Config {

   /**
    * @brief Implementation of the XML configuration file reader
    */
   class ConfigurationReader: public IConfigurationFile<Io::Xml::IXmlReader>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim  Dimensionality of simulation
          * @param type Type of the simulation
          */
         ConfigurationReader(const int dim, const std::vector<bool>& isPeriodicBox, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ConfigurationReader();

         /**
          * @brief Read content of configuration file
          */
         virtual void read();
         
      protected:
         /**
          * @brief Read configuration block
          *
          * @param required Configuration tags are required?
          */
         void readBlock(SharedIConfigurationBlock spBlock, const bool required);

      private:
         /**
          * @brief NAME of the configuration file
          */
         static const std::string NAME;
   };

   /// Typedef for a smart pointer of a ConfigurationReader
   typedef std::shared_ptr<ConfigurationReader> SharedConfigurationReader;

}
}
}

#endif // QUICC_IO_CONFIG_CONFIGURATIONREADER_HPP
