/** 
 * @file ShellTorPolLSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics L energy spectrum calculation for a Toroidal/Poloidal field in a spherical shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLTORPOLLSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLTORPOLLSPECTRUMWRITER_HPP

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
#include "QuICC/Io/Variable/ISphericalTorPolLSpectrumWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics L energy spectrum calculation for a Toroidal/Poloidal field in a spherical shell
    */
   class ShellTorPolLSpectrumWriter: public ISphericalTorPolLSpectrumWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ShellTorPolLSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ShellTorPolLSpectrumWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<ShellTorPolLSpectrumWriter> SharedShellTorPolLSpectrumWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_SHELLTORPOLLSPECTRUMWRITER_HPP
