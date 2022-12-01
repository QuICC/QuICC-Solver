/** 
 * @file ISphericalTorPolRSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics radial power spectrum calculation for a Toroidal/Poloidal field in a spherical geometry
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALTORPOLRSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALTORPOLRSPECTRUMWRITER_HPP

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
    * @brief Implementation of the ASCII spherical harmonics radial power spectrum calculation for a Toroidal/Poloidal field in a spherical geometry
    */
   class ISphericalTorPolRSpectrumWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          */
         ISphericalTorPolRSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalTorPolRSpectrumWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init() override;

         /**
          * @brief Activate output of parity splitting in power output
          */
         void showParity();

         /**
          * @brief Compute power
          */
         void compute(Transform::TransformCoordinatorType& coord) override;

         /**
          * @brief Requires heavy calculation?
          */
         virtual bool isHeavy() const override;
         
      protected:
         /**
          * @brief Prepare spectral field data for computation
          */
         void prepareInput(const FieldComponents::Spectral::Id sId, Transform::TransformCoordinatorType& coord);

         /**
          * @brief Write content
          */
         virtual void writeContent() override;

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
          * @brief Radial grid
          */
         Array mGrid;

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
         virtual void resetPower();

         /**
          * @brief Store power from Q component
          */
         virtual void storeQPower(const int n, const int l, const int m, const MHDFloat power);

         /**
          * @brief Store power from S component
          */
         virtual void storeSPower(const int n, const int l, const int m, const MHDFloat power);

         /**
          * @brief Store power from T component
          */
         virtual void storeTPower(const int n, const int l, const int m, const MHDFloat power);
   };

   inline bool ISphericalTorPolRSpectrumWriter::isHeavy() const
   {
      return true;
   }

}
}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALTORPOLPOWERBASEWRITER_HPP
