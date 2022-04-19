/** 
 * @file SphereTorPolEnstrophyLSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics enstrophy L Spectrum calculation for a Toroidal/Poloidal field in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERETORPOLENSTROPHYLSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERETORPOLENSTROPHYLSPECTRUMWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Io/Variable/ISphericalTorPolEnstrophyLSpectrumWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics enstrophy L spectrum calculation for a Toroidal/Poloidal field in a sphere
    */
   class SphereTorPolEnstrophyLSpectrumWriter: public ISphericalTorPolEnstrophyLSpectrumWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         SphereTorPolEnstrophyLSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~SphereTorPolEnstrophyLSpectrumWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer
   typedef std::shared_ptr<SphereTorPolEnstrophyLSpectrumWriter> SharedSphereTorPolEnstrophyLSpectrumWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_SPHERETORPOLENSTROPHYLSPECTRUMWRITER_HPP
