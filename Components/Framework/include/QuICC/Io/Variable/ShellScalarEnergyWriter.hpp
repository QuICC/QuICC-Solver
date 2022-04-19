/** 
 * @file ShellScalarEnergyWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a scalar field in a spherical shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLSCALARENERGYWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLSCALARENERGYWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Io/Variable/ISphericalScalarEnergyWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics energy calculation for a scalar field in a spherical shell
    */
   class ShellScalarEnergyWriter: public ISphericalScalarEnergyWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ShellScalarEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ShellScalarEnergyWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<ShellScalarEnergyWriter> SharedShellScalarEnergyWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_SHELLSCALARENERGYWRITER_HPP
