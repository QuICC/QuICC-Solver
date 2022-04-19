/** 
 * @file ShellScalarLSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics L energy spectrum calculation for a scalar field in a spherical shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLSCALARLSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLSCALARLSPECTRUMWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Io/Variable/ISphericalScalarLSpectrumWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics L energy spectrum calculation for a scalar field in a spherical shell
    */
   class ShellScalarLSpectrumWriter: public ISphericalScalarLSpectrumWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ShellScalarLSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ShellScalarLSpectrumWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<ShellScalarLSpectrumWriter> SharedShellScalarLSpectrumWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_SHELLSCALARLSPECTRUMWRITER_HPP
