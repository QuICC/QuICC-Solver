/** 
 * @file ShellTorPolEnstrophyWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics enstrophy calculation for a Toroidal/Poloidal field in a spherical shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLTORPOLENSTROPHYWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLTORPOLENSTROPHYWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Io/Variable/ISphericalTorPolEnstrophyWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics enstrophy calculation for a Toroidal/Poloidal field in a spherical shell
    */
   class ShellTorPolEnstrophyWriter: public ISphericalTorPolEnstrophyWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ShellTorPolEnstrophyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ShellTorPolEnstrophyWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<ShellTorPolEnstrophyWriter> SharedShellTorPolEnstrophyWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_SHELLTORPOLENSTROPHYWRITER_HPP
