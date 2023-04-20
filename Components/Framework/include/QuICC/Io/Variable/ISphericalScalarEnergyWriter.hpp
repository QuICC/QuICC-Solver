/** 
 * @file ISphericalScalarEnergyWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a scalar field in a spherical geometry
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALSCALARENERGYWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALSCALARENERGYWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Io/Variable/ISphericalScalarEnergyBaseWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics energy calculation for a scalar field in a spherical geometry
    */
   class ISphericalScalarEnergyWriter: public ISphericalScalarEnergyBaseWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ISphericalScalarEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalScalarEnergyWriter() = default;
         
      protected:
         /**
          * @brief Write content
          */
         virtual void writeContent();

      private:
         /**
          * @brief Storage for the scalar energy
          */
         Array mEnergy;

         /**
          * @brief Reset energy storage
          */
         virtual void resetEnergy();

         /**
          * @brief Store energy
          *
          * @param l      Harmonic degree
          * @param m      Harmonic order
          * @param energy  Energy of mode
          */
         virtual void storeEnergy(const int l, const int m, const MHDFloat energy);
   };

} // Variable
} // Io
} // QuICC

#endif // QUICC_IO_VARIABLE_ISPHERICALSCALARENERGYWRITER_HPP
