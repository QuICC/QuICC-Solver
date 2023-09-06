/** 
 * @file ISphericalTorPolEnergyWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a spherical geometry
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALTORPOLENERGYWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALTORPOLENERGYWRITER_HPP

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
    * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a spherical geometry
    */
   class ISphericalTorPolEnergyWriter: public ISphericalTorPolEnergyBaseWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ISphericalTorPolEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalTorPolEnergyWriter() = default;
         
      protected:
         /**
          * @brief Write content
          */
         virtual void writeContent();

      private:
         /**
          * @brief Storage for the Toroidal energy
          */
         Array mTorEnergy;

         /**
          * @brief Storage for the Poloidal energy
          */
         Array mPolEnergy;

         /**
          * @brief Reset energy storage
          */
         virtual void resetEnergy();

         /**
          * @brief Store energy from Q component
          *
          * @param l      Harmonic degree
          * @param m      Harmonic order
          * @param energy  Energy of mode
          */
         virtual void storeQEnergy(const int l, const int m, const MHDFloat energy);

         /**
          * @brief Store energy from S component
          *
          * @param l      Harmonic degree
          * @param m      Harmonic order
          * @param energy  Energy of mode
          */
         virtual void storeSEnergy(const int l, const int m, const MHDFloat energy);

         /**
          * @brief Store energy from T component
          *
          * @param l      Harmonic degree
          * @param m      Harmonic order
          * @param energy  Energy of mode
          */
         virtual void storeTEnergy(const int l, const int m, const MHDFloat energy);
   };

} // Variable
} // Io
} // QuICC

#endif // QUICC_IO_VARIABLE_ISPHERICALTORPOLENERGYWRITER_HPP
