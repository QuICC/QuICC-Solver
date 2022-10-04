/** 
 * @file ISphericalTorPolEnstrophyBaseWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics enstrophy calculation for a Toroidal/Poloidal field in a spherical geometry
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALTORPOLENSTROPHYWBASERITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALTORPOLENSTROPHYWBASERITER_HPP

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
    * @brief Implementation of the ASCII spherical harmonics enstrophy calculation for a Toroidal/Poloidal field in a spherical geometry
    */
   class ISphericalTorPolEnstrophyBaseWriter: public IVariableAsciiWriter
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
         ISphericalTorPolEnstrophyBaseWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const IAsciiWriter::WriteMode mode = IAsciiWriter::EXTEND);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalTorPolEnstrophyBaseWriter();

         /**
          * @brief Activate output of parity splitting in enstrophy output
          */
         void showParity();

         /**
          * @brief Compute enstrophy
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
          * @brief Spherical volume to normalize enstrophy to enstrophy density
          */
         MHDFloat mVolume;

         /**
          * @brief Flag to show parity split in enstrophy
          */
         bool mShowParity;

      private:
         /**
          * @brief Reset enstrophy storage
          */
         virtual void resetEnstrophy() = 0;

         /**
          * @brief Store enstrophy from Q component
          */
         virtual void storeQEnstrophy(const int l, const int m, const MHDFloat enstrophy) = 0;

         /**
          * @brief Store enstrophy from S component
          */
         virtual void storeSEnstrophy(const int l, const int m, const MHDFloat enstrophy) = 0;

         /**
          * @brief Store enstrophy from T component
          */
         virtual void storeTEnstrophy(const int l, const int m, const MHDFloat enstrophy) = 0;
   };

   inline bool ISphericalTorPolEnstrophyBaseWriter::isHeavy() const
   {
      return true;
   }

}
}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALTORPOLENSTROPHYBASEWRITER_HPP
