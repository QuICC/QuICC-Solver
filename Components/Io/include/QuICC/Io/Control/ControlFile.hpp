/** 
 * @file ControlFile.hpp
 * @brief Implementation of a control file
 */

#ifndef QUICC_IO_CONTROL_CONTROLFILE_HPP
#define QUICC_IO_CONTROL_CONTROLFILE_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Io {

namespace Control {

   /**
    * @brief Implementation of a control file
    */
   class ControlFile
   {
      public:
         /**
         * @brief Constructor
         */
         ControlFile();

         /**
         * @brief Destructor
         */
         virtual ~ControlFile();

         /**
          * @brief Get filename
          */
         std::string  filename() const;
         
      protected:
         /**
          * @brief Get the full header
          */
         std::string  header() const;

         /**
          * @brief Get the full type
          */
         std::string  type() const;

         /**
          * @brief Get the full version
          */
         std::string  version() const;

         /**
          * @brief Get the name
          */
         std::string  name() const;

         /**
          * @brief Get extension
          */
         std::string  extension() const;

      private:
         /**
          * @brief Header of file before header information
          */
         static const std::string FILE_HEADER;

         /**
          * @brief Header of file before optional information
          */
         static const std::string FILE_TYPE;

         /**
          * @brief Header of file before version information
          */
         static const std::string FILE_VERSION;

         /**
          * @brief NAME of the runtime control file
          */
         static const std::string FILE_NAME;

         /**
          * @brief EXTENSION of the runtime control file
          */
         static const std::string FILE_EXTENSION;
   };

}
}
}

#endif // QUICC_IO_CONTROL_CONTROLFILE_HPP
