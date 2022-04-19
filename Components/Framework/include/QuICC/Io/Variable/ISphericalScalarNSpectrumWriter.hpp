/** 
 * @file ISphericalScalarNSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics L power spectrum calculation for a scalar field in a spherical geometry
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALSCALARNSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALSCALARNSPECTRUMWRITER_HPP

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
#include "QuICC/Io/Variable/ISphericalScalarPowerBaseWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics L power spectrum calculation for a scalar field in a spherical geometry
    */
   class ISphericalScalarNSpectrumWriter: public ISphericalScalarPowerBaseWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ISphericalScalarNSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalScalarNSpectrumWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init() override;
         
      protected:
         /**
          * @brief Write content
          */
         virtual void writeContent() override;

      private:
         /**
          * @brief Storage for the scalar power
          */
         Matrix mPower;

         /**
          * @brief Initialize power storage
          */
         virtual void initializePower() override;

         /**
          * @brief Store power
          */
         virtual void storePower(const int n, const int l, const int m, const MHDFloat power) override;
   };

}
}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALSCALARNSPECTRUMWRITER_HPP
