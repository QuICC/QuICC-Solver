/**
 * @file DirectAsciiWriter.hpp
 * @brief Implementation of a direct access ASCII writer
 */

#ifndef QUICC_IO_ASCII_DIRECTASCIIWRITER_HPP
#define QUICC_IO_ASCII_DIRECTASCIIWRITER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Io/Ascii/IAsciiWriter.hpp"

namespace QuICC {

namespace Io {

namespace Ascii {

   /**
    * @brief Implementation of a direct access ASCII writer
    */
   class DirectAsciiWriter: public IAsciiWriter
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
         DirectAsciiWriter(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~DirectAsciiWriter();

         /**
          * @brief Initialise the file
          */
         virtual void init();

         /**
          * @brief Initialise the file in debug mode (no IO filter)
          */
         virtual void initDebug();

         /**
          * @brief Initialise the file without header
          */
         void initNoHeader();

         /**
          * @brief Finalise the file
          */
         virtual void finalize();

         /**
          * @brief Finalise the file in debug mode (no IO filter)
          */
         virtual void finalizeDebug();

         /**
          * @brief Get direct access to file handle
          */
         std::ofstream& file();

      protected:
         /**
          * @brief This call does nothing in this case
          */
         virtual void writeContent() {};

      private:
   };
}
}
}

#endif // QUICC_IO_ASCII_DIRECTASCIIWRITER_HPP
