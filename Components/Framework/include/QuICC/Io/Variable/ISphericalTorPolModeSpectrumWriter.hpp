/** 
 * @file ISphericalTorPolModeSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics mode energy spectrum calculation for a Toroidal/Poloidal field in a spherical geometry
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALTORPOLMODESPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALTORPOLMODESPECTRUMWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Io/Variable/ISphericalTorPolEnergyBaseWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics mode energy spectrum calculation for a Toroidal/Poloidal field in a spherical geometry
    */
   class ISphericalTorPolModeSpectrumWriter: public ISphericalTorPolEnergyBaseWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ISphericalTorPolModeSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalTorPolModeSpectrumWriter() = default;

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:
         /**
          * @brief Write content
          */
         virtual void writeContent();

      private:
         /**
          * @brief Storage for the Toroidal energy
          */
         Matrix mTorEnergy;

         /**
          * @brief Storage for the Poloidal energy
          */
         Matrix mPolEnergy;

         /**
          * @brief Reset energy storage
          */
         virtual void resetEnergy();

         /**
          * @brief Store energy from Q component
          *
          * @param l       harmonic degree
          * @param m       harmonic order
          * @param energy  energy in (l,m) mode
          */
         virtual void storeQEnergy(const int l, const int m, const MHDFloat energy);

         /**
          * @brief Store energy from S component
          *
          * @param l       harmonic degree
          * @param m       harmonic order
          * @param energy  energy in (l,m) mode
          */
         virtual void storeSEnergy(const int l, const int m, const MHDFloat energy);

         /**
          * @brief Store energy from T component
          *
          * @param l       harmonic degree
          * @param m       harmonic order
          * @param energy  energy in (l,m) mode
          */
         virtual void storeTEnergy(const int l, const int m, const MHDFloat energy);
   };

}
}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALTORPOLMODESPECTRUMWRITER_HPP
