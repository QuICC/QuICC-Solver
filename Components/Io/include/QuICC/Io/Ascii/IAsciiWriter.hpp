/**
 * @file IAsciiWriter.hpp
 * @brief General interface to an ASCII writer
 */

#ifndef QUICC_IO_ASCII_IASCIIWRITER_HPP
#define QUICC_IO_ASCII_IASCIIWRITER_HPP

// Configuration includes
//

// System includes
//
#include <fstream>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Io/Ascii/AsciiFile.hpp"

namespace QuICC {

namespace Io {

namespace Ascii {

   /**
    * @brief General interface to an ASCII writer
    */
   class IAsciiWriter: public AsciiFile
   {
      public:
         /**
          * @brief Possible write modes
          */
         enum WriteMode {
            EXTEND, // New values are appended
            OVERWRITE, // File is overwritten
            NUMBER, // New file is created with incremented number
         };

         /**
         * @brief Constructor
         *
         * @param name     Filename
         * @param ext      File extension
         * @param header   Header string of file
         * @param type     Type string of file
         * @param version  Version string of file
         * @param mode     Write mode of file
         */
         IAsciiWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const WriteMode mode = EXTEND);

         /**
         * @brief Destructor
         */
         virtual ~IAsciiWriter();

         /**
          * @brief Initialise the file
          */
         virtual void init();

         /**
          * @brief Write the content
          */
         void write();

         /**
          * @brief Finalise the file
          */
         virtual void finalize();

         /**
          * @brief Set output mode to OVERWRITE the file
          */
         void overwriteOutput();

         /**
          * @brief Set output mode to create new numbered file
          */
         void numberOutput();

         /**
          * @brief Set output to extend data in file
          */
         void extendOutput();

         /**
          * @brief Set output frequency of numbered file
          */
         void onlyEvery(const int n);

         /**
          * @brief Is file active
          */
         bool isActive() const;

      protected:
         /**
          * @brief Operation to perform just before writing data
          */
         void preWrite();

         /**
          * @brief Write implementation called by write
          */
         virtual void writeContent() = 0;

         /**
          * @brief Operation to perform just after writing data
          */
         void postWrite();

         /**
          * @brief Handle to the file
          */
         std::ofstream mFile;

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Open the file in debug mode (no IO filter)
          */
         void openDebug();

         /**
          * @brief Close the file
          */
         void close();

         /**
          * @brief Close the file in debug mode (no IO filter)
          */
         void closeDebug();

         /**
          * @brief Create the file
          */
         void create();

         /**
          * @brief end the file
          */
         void end();

         /**
          * @brief Update the file name with counter value
          */
         void updateName();

      private:
         /**
          * @brief Zero Fill width
          */
         static const int msIDWidth;

         /**
          * @brief Is file initialized?
          */
         bool mIsInitialized;

         /**
          * @brief File counter
          */
         int mCounter;

         /**
          * @brief Only output every X file
          */
         int mEvery;

         /**
          * @brief Base name of the file
          */
         const std::string    mBaseName;

         /**
          * @brief Write mode of file
          */
         WriteMode mMode;
   };

   /// Typedef for a shared pointer of a IAsciiWriter
   typedef std::shared_ptr<IAsciiWriter> SharedIAsciiWriter;
}
}
}

#endif // QUICC_IO_ASCII_IASCIIWRITER_HPP
