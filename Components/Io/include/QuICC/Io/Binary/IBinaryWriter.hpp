/** 
 * @file IBinaryWriter.hpp
 * @brief Interface to a general binary writer
 */

#ifndef QUICC_IO_BINARY_IBINARYWRITER_HPP
#define QUICC_IO_BINARY_IBINARYWRITER_HPP

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
#include "QuICC/Io/Binary/BinaryFile.hpp"

namespace QuICC {

namespace Io {

namespace Binary {

   /**
    * @brief Interface to a general binary writer
    */
   class IBinaryWriter: public BinaryFile
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
         IBinaryWriter(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~IBinaryWriter();

         /**
          * @brief Initialise the file
          */
         virtual void init() = 0;

         /**
          * @brief Write the content
          */
         virtual void write() = 0;

         /**
          * @brief Finalise the file
          */
         virtual void finalize() = 0;
         
      protected:

         /**
          * @brief Handle to the file
          */
         std::ofstream mFile;

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Close the file
          */
         void close();

      private:
   };

   /// Typedef for a shared pointer of a IBinaryWriter
   typedef std::shared_ptr<IBinaryWriter> SharedIBinaryWriter;
}
}
}

#endif // QUICC_IO_BINARY_IBINARYWRITER_HPP
