/** 
 * @file ISphericalTorPolPowerBaseWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics power spectrum calculation for a Toroidal/Poloidal field in a spherical geometry
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALTORPOLPOWERBASEWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALTORPOLPOWERBASEWRITER_HPP

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
#include "QuICC/Io/Variable/IVariableAsciiWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics power spectrum calculation for a Toroidal/Poloidal field in a spherical geometry
    */
   class ISphericalTorPolPowerBaseWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name       Filename
          * @param ext        File extension
          * @param header     Header string of file
          * @param type       Type string of file
          * @param version    Version string of file
          * @param id         ID of the dimension space
          * @param mode       Write mode of file
          */
         ISphericalTorPolPowerBaseWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const IAsciiWriter::WriteMode mode = IAsciiWriter::EXTEND);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalTorPolPowerBaseWriter();

         /**
          * @brief Activate output of parity splitting in power output
          */
         void showParity();

         /**
          * @brief Compute power
          */
         void compute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief Requires heavy calculation?
          */
         virtual bool isHeavy() const;
         
      protected:
         /**
          * @brief Prepare spectral field data for computation
          */
         void prepareInput(const FieldComponents::Spectral::Id sId, Transform::TransformCoordinatorType& coord);

         /**
          * @brief Data ordering is m slowest
          */
         bool mHasMOrdering;

         /**
          * @brief Spherical volume to normalize power to power density
          */
         MHDFloat mVolume;

         /**
          * @brief Flag to show parity split in power
          */
         bool mShowParity;

      private:
         /**
          * @brief Reset power storage
          */
         virtual void resetPower() = 0;

         /**
          * @brief Store power from Q component
          */
         virtual void storeQPower(const int n, const int l, const int m, const MHDFloat power) = 0;

         /**
          * @brief Store power from S component
          */
         virtual void storeSPower(const int n, const int l, const int m, const MHDFloat power) = 0;

         /**
          * @brief Store power from T component
          */
         virtual void storeTPower(const int n, const int l, const int m, const MHDFloat power) = 0;
   };

   inline bool ISphericalTorPolPowerBaseWriter::isHeavy() const
   {
      return true;
   }

}
}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALTORPOLPOWERBASEWRITER_HPP
