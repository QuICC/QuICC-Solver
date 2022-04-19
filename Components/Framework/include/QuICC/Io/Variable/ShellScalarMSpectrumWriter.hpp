/** 
 * @file ShellScalarMSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics M energy spectrum calculation for a scalar field in a spherical shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLSCALARMSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLSCALARMSPECTRUMWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Io/Variable/ISphericalScalarMSpectrumWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics M energy spectrum calculation for a scalar field in a spherical shell
    */
   class ShellScalarMSpectrumWriter: public ISphericalScalarMSpectrumWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ShellScalarMSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ShellScalarMSpectrumWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<ShellScalarMSpectrumWriter> SharedShellScalarMSpectrumWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_SHELLSCALARMSPECTRUMWRITER_HPP
