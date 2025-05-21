/** 
 * @file ShellTorPolModeSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics mode energy spectrum calculation for a Toroidal/Poloidal field in a spherical shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLTORPOLMODESPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLTORPOLMODESPECTRUMWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Io/Variable/ISphericalTorPolModeSpectrumWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics mode energy spectrum calculation for a Toroidal/Poloidal field in a spherical shell
    */
   class ShellTorPolModeSpectrumWriter: public ISphericalTorPolModeSpectrumWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ShellTorPolModeSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ShellTorPolModeSpectrumWriter() = default;

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<ShellTorPolModeSpectrumWriter> SharedShellTorPolModeSpectrumWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_SHELLTORPOLMODESPECTRUMWRITER_HPP
