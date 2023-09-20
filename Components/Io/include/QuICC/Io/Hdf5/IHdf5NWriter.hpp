/**
 * @file IHdf5NWriter.hpp
 * @brief Interface to a numbering HDF5 file
 */

#ifndef QUICC_IO_HDF5_IHDF5NWRITER_HPP
#define QUICC_IO_HDF5_IHDF5NWRITER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Io/Hdf5/IHdf5Writer.hpp"

namespace QuICC {

namespace Io {

namespace Hdf5 {

   /**
    * @brief Class describes an hdf5 file writer where the data gets written to a new file every time
    *          a number is appended to the name.
    */
   class IHdf5NWriter: public IHdf5Writer
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
         IHdf5NWriter(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
          * @brief Destructor
          */
         virtual ~IHdf5NWriter();

         /**
          * @brief Initialise the file
          */
         void init();

         /**
          * @brief Write the content
          */
         virtual void write() = 0;

         /**
          * @brief Finalise the file
          */
         virtual void finalize();

         /**
          * @brief Change the basename
          *
          * @param base New basename
          */
         void changeBasename(std::string base);

      protected:
         /**
          * @brief Update the filename with new number
          */
         void updateName();

         /**
          * @brief Operation to perform just before writing
          */
         void preWrite();

         /**
          * @brief Operations to perform just after writing
          */
         void postWrite();

      private:
         /**
          * @brief Width of the zero fill for the file id
          */
         static const int ZEROWIDTH;

         /**
          * @brief Counter for the file number
          */
         int mCounter;

         /**
          * @brief Base name used for appending the number
          */
         std::string mBaseName;
   };

}
}
}

#endif // QUICC_IO_HDF5_IHDF5WRITER_HPP
