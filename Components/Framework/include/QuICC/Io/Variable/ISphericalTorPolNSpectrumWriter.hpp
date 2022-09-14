/** 
 * @file ISphericalTorPolNSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics power calculation for a Toroidal/Poloidal field in a spherical geometry
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALTORPOLNSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALTORPOLNSPECTRUMWRITER_HPP

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
#include "QuICC/Io/Variable/ISphericalTorPolPowerBaseWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics power calculation for a Toroidal/Poloidal field in a spherical geometry
    */
   class ISphericalTorPolNSpectrumWriter: public ISphericalTorPolPowerBaseWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ISphericalTorPolNSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalTorPolNSpectrumWriter();

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
          * @brief Storage for the Toroidal power
          */
         Matrix mTorPower;

         /**
          * @brief Storage for the Poloidal power
          */
         Matrix mPolPower;

         /**
          * @brief Reset power storage
          */
         virtual void resetPower() override;

         /**
          * @brief Store power from Q component
          */
         virtual void storeQPower(const int n, const int l, const int m, const MHDFloat power) override;

         /**
          * @brief Store power from S component
          */
         virtual void storeSPower(const int n, const int l, const int m, const MHDFloat power) override;

         /**
          * @brief Store power from T component
          */
         virtual void storeTPower(const int n, const int l, const int m, const MHDFloat power) override;
   };

}
}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALTORPOLNSPECTRUMWRITER_HPP
