/** 
 * @file SphereTorPolEnstrophyWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics enstrophy calculation for a Toroidal/Poloidal field in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERETORPOLENSTROPHYWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERETORPOLENSTROPHYWRITER_HPP

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
    * @brief Implementation of the ASCII spherical harmonics enstrophy calculation for a Toroidal/Poloidal field in a sphere
    */
   class SphereTorPolEnstrophyWriter: public ISphericalTorPolEnstrophyWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         SphereTorPolEnstrophyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~SphereTorPolEnstrophyWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<SphereTorPolEnstrophyWriter> SharedSphereTorPolEnstrophyWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_SPHERETORPOLENSTROPHYWRITER_HPP
