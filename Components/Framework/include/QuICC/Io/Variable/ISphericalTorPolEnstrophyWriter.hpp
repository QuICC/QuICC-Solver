/** 
 * @file ISphericalTorPolEnstrophyWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics enstrophy calculation for a Toroidal/Poloidal field in a spherical geometry
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALTORPOLENSTROPHYWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALTORPOLENSTROPHYWRITER_HPP

// Configuration includes
//
#include <memory>

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Io/Variable/ISphericalTorPolEnstrophyBaseWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics enstrophy calculation for a Toroidal/Poloidal field in a spherical geometry
    */
   class ISphericalTorPolEnstrophyWriter: public ISphericalTorPolEnstrophyBaseWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ISphericalTorPolEnstrophyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalTorPolEnstrophyWriter();
         
      protected:
         /**
          * @brief Write content
          */
         virtual void writeContent();

      private:
         /**
          * @brief Storage for the Toroidal enstrophy
          */
         Array mTorEnstrophy;

         /**
          * @brief Storage for the Poloidal enstrophy
          */
         Array mPolEnstrophy;

         /**
          * @brief Reset enstrophy storage
          */
         virtual void resetEnstrophy();

         /**
          * @brief Store enstrophy from Q component
          */
         virtual void storeQEnstrophy(const int l, const int m, const MHDFloat enstrophy);

         /**
          * @brief Store enstrophy from S component
          */
         virtual void storeSEnstrophy(const int l, const int m, const MHDFloat enstrophy);

         /**
          * @brief Store enstrophy from T component
          */
         virtual void storeTEnstrophy(const int l, const int m, const MHDFloat enstrophy);
   };

}
}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALTORPOLENSTROPHYWRITER_HPP
