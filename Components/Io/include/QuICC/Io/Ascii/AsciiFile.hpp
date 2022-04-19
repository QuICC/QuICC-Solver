/**
 * @file AsciiFile.hpp
 * @brief Implementation of a general ASCII file
 */

#ifndef QUICC_IO_ASCII_ASCIIFILE_HPP
#define QUICC_IO_ASCII_ASCIIFILE_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Io {

namespace Ascii {

   /**
    * @brief Implementation of a general ASCII file
    */
   class AsciiFile
   {
      public:
         /**
         * @brief Constructor
         *
         * @param name     Filename
         * @param ext      File extension
         * @param header   Header string of file
         * @param type     Type string of file
         * @param version  Version string of file
         */
         AsciiFile(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~AsciiFile();

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

         /**
          * @brief Reset name
          *
          * @param name New name
          */
         void resetName(std::string name);

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
          * @brief Name of the file without extension
          */
         std::string mName;

         /**
          * @brief File extension
          */
         std::string mExt;

         /**
          * @brief Header of the file to check compatibility
          */
         std::string mHeader;

         /**
          * @brief Type of the file to check compatibility
          */
         std::string mType;

         /**
          * @brief Version of the file to check compatibility
          */
         std::string mVersion;
   };

}
}
}

#endif // QUICC_IO_ASCII_ASCIIFILE_HPP
