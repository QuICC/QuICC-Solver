/** 
 * @file ISphericalScalarLSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics L energy spectrum calculation for a scalar field in a spherical geometry
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALSCALARLSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALSCALARLSPECTRUMWRITER_HPP

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
#include "QuICC/Io/Variable/ISphericalScalarEnergyBaseWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics L energy spectrum calculation for a scalar field in a spherical geometry
    */
   class ISphericalScalarLSpectrumWriter: public ISphericalScalarEnergyBaseWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ISphericalScalarLSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalScalarLSpectrumWriter();

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
          * @brief Storage for the scalar energy
          */
         Array mEnergy;

         /**
          * @brief Initialize energy storage
          */
         virtual void initializeEnergy();

         /**
          * @brief Store energy
          */
         virtual void storeEnergy(const int l, const int m, const MHDFloat energy);
   };

}
}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALSCALARLSPECTRUMWRITER_HPP
