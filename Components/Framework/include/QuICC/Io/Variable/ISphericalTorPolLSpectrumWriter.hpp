/** 
 * @file ISphericalTorPolLSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a spherical geometry
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALTORPOLLSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALTORPOLLSPECTRUMWRITER_HPP

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
#include "QuICC/Io/Variable/ISphericalTorPolEnergyBaseWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a spherical geometry
    */
   class ISphericalTorPolLSpectrumWriter: public ISphericalTorPolEnergyBaseWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ISphericalTorPolLSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalTorPolLSpectrumWriter();

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
         Array mTorEnergy;

         /**
          * @brief Storage for the Poloidal energy
          */
         Array mPolEnergy;

         /**
          * @brief Initialize energy storage
          */
         virtual void initializeEnergy();

         /**
          * @brief Store energy from Q component
          */
         virtual void storeQEnergy(const int l, const int m, const MHDFloat energy);

         /**
          * @brief Store energy from S component
          */
         virtual void storeSEnergy(const int l, const int m, const MHDFloat energy);

         /**
          * @brief Store energy from T component
          */
         virtual void storeTEnergy(const int l, const int m, const MHDFloat energy);
   };

}
}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALTORPOLLSPECTRUMWRITER_HPP
