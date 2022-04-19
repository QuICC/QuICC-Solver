/**
 * @file IAsciiReader.hpp
 * @brief Interface to an ASCII file reader
 */

#ifndef QUICC_IO_ASCII_IASCIIREADER_HPP
#define QUICC_IO_ASCII_IASCIIREADER_HPP

// System includes
//
#include <fstream>

// External includes
//

// Project includes
//
#include "QuICC/Io/Ascii/AsciiFile.hpp"

namespace QuICC {

namespace Io {

namespace Ascii {

   /**
    * @brief Interface to an ASCII file reader
    */
   class IAsciiReader: public AsciiFile
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
         IAsciiReader(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~IAsciiReader();

         /**
          * @brief Initialise the file
          */
         virtual void init();

         /**
          * @brief Read the content
          */

         virtual void read() = 0;

         /**
          * @brief Finalise the file
          */
         virtual void finalize();

      protected:
         /**
          * @brief Handle to the file
          */
         std::ifstream mFile;

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Close the file
          */
         void close();

         /**
          * @brief Check compatibility of opened file
          */
         virtual void checkCompatibility();

      private:
   };

}
}
}

#endif // QUICC_IO_ASCII_IASCIIREADER_HPP
