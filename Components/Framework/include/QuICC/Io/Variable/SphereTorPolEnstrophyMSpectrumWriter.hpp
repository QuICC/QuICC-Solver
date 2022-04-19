/** 
 * @file SphereTorPolEnstrophyMSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics enstrophy M spectrum calculation for a Toroidal/Poloidal field in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERETORPOLENSTROPHYMSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERETORPOLENSTROPHYMSPECTRUMWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Io/Variable/ISphericalTorPolEnstrophyMSpectrumWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics enstrophy M spectrum calculation for a Toroidal/Poloidal field in a sphere
    */
   class SphereTorPolEnstrophyMSpectrumWriter: public ISphericalTorPolEnstrophyMSpectrumWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         SphereTorPolEnstrophyMSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~SphereTorPolEnstrophyMSpectrumWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer
   typedef std::shared_ptr<SphereTorPolEnstrophyMSpectrumWriter> SharedSphereTorPolEnstrophyMSpectrumWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_SPHERETORPOLENSTROPHYMSPECTRUMWRITER_HPP
