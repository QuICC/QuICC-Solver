/** 
 * @file ISphericalTorPolEnergyBaseWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a spherical geometry
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALTORPOLENERGYWBASERITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALTORPOLENERGYWBASERITER_HPP

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
    * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a spherical geometry
    */
   class ISphericalTorPolEnergyBaseWriter: public IVariableAsciiWriter
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
         ISphericalTorPolEnergyBaseWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const IAsciiWriter::WriteMode mode = IAsciiWriter::EXTEND);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalTorPolEnergyBaseWriter();

         /**
          * @brief Activate output of parity splitting in energy output
          */
         void showParity();

         /**
          * @brief Compute energy
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
          * @brief Spherical volume to normalize energy to energy density
          */
         MHDFloat mVolume;

         /**
          * @brief Flag to show parity split in energy
          */
         bool mShowParity;

      private:
         /**
          * @brief Reset energy storage
          */
         virtual void resetEnergy() = 0;

         /**
          * @brief Store energy from Q component
          */
         virtual void storeQEnergy(const int l, const int m, const MHDFloat energy) = 0;

         /**
          * @brief Store energy from S component
          */
         virtual void storeSEnergy(const int l, const int m, const MHDFloat energy) = 0;

         /**
          * @brief Store energy from T component
          */
         virtual void storeTEnergy(const int l, const int m, const MHDFloat energy) = 0;
   };

   inline bool ISphericalTorPolEnergyBaseWriter::isHeavy() const
   {
      return true;
   }

}
}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALTORPOLENERGYBASEWRITER_HPP
