/** 
 * @file ShellTorPolEnergyWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a spherical shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLTORPOLENERGYWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLTORPOLENERGYWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Io/Variable/ISphericalTorPolEnergyWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a spherical shell
    */
   class ShellTorPolEnergyWriter: public ISphericalTorPolEnergyWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ShellTorPolEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ShellTorPolEnergyWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<ShellTorPolEnergyWriter> SharedShellTorPolEnergyWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_SHELLTORPOLENERGYWRITER_HPP
