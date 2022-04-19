/**
 * @file StdOutPipe.hpp
 * @brief Implementation of standard output pipe file
 */

#ifndef QUICC_IO_ASCII_STDOUTPIPE_HPP
#define QUICC_IO_ASCII_STDOUTPIPE_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Io/Ascii/IAsciiWriter.hpp"

namespace QuICC {

namespace Io {

namespace Ascii {

   /**
    * @brief Implementation of the standard output pipe into an ASCII file
    */
   class StdOutPipe: public IAsciiWriter
   {
      public:
         /**
         * @brief Constructor
         *
         * @param name Name of the std message ouput file
         */
         StdOutPipe(std::string name);

         /**
         * @brief Destructor
         */
         virtual ~StdOutPipe();

         /**
          * @brief Init the file
          */
         void init();

         /**
          * @brief Finalise
          */
         void finalize();

      protected:
         /**
          * @brief
          */
         virtual void writeContent();

      private:
         /**
          * @brief HEADER part for StdOutPipe file
          */
         static const std::string  HEADER;

         /**
          * @brief TYPE part for StdOutPipe file
          */
         static const std::string  TYPE;

         /**
          * @brief VERSION part for StdOutPipe file
          */
         static const std::string  VERSION;

         /**
          * @brief BASENAME of StdOutPipe file
          */
         static const std::string  BASENAME;

         /**
          * @brief EXTENSION of StdOutPipe file
          */
         static const std::string  EXTENSION;

         /**
          * Backup std::cout buffer
          */
         std::streambuf   *mpCoutBuffer;
   };

   /// Typedef for a smart shared pointer of a StdOutPipe
   typedef std::shared_ptr<StdOutPipe> SharedStdOutPipe;

}
}
}

#endif // QUICC_IO_ASCII_STDOUTPIPE_HPP
