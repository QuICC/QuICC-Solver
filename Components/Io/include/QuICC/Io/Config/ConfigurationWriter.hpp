/**
 * @file ConfigurationWriter.hpp
 * @brief Implementation of the XML configuration file writer
 */

#ifndef QUICC_IO_CONFIG_CONFIGURATIONWRITER_HPP
#define QUICC_IO_CONFIG_CONFIGURATIONWRITER_HPP

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
#include "QuICC/Io/Xml/IXmlWriter.hpp"
#include "QuICC/Io/Config/IConfigurationFile.hpp"

namespace QuICC {

namespace Io {

namespace Config {

   /**
    * @brief Implementation of the XML configuration file writer
    */
   class ConfigurationWriter: public IConfigurationFile<Io::Xml::IXmlWriter>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim  Dimensionality of simulation
          * @param type Type of the simulation
          */
         ConfigurationWriter(const int dim, const std::vector<bool>& isPeriodicBox, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ConfigurationWriter();

         /**
          * @brief Read content of configuration file
          */
         virtual void write();

      protected:
         /**
          * @brief Write configuration block
          */
         void writeBlock(SharedIConfigurationBlock spBlock);

         /**
          * @brief Write XML root
          */
         void writeRoot();

      private:
         /**
          * @brief NAME of the configuration file
          */
         static const std::string NAME;
   };

   /// Typedef for a smart pointer of a ConfigurationWriter
   typedef std::shared_ptr<ConfigurationWriter> SharedConfigurationWriter;

}
}
}

#endif // QUICC_IO_CONFIG_CONFIGURATIONWRITER_HPP
