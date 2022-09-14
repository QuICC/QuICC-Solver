/** 
 * @file ISphericalScalarRSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics radial power spectrum calculation for a scalar field in a spherical geometry
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALSCALARRARSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALSCALARRARSPECTRUMWRITER_HPP

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
    * @brief Implementation of the ASCII spherical harmonics radial power spectrum calculation for a scalar field in a spherical geometry
    */
   class ISphericalScalarRSpectrumWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ISphericalScalarRSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalScalarRSpectrumWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init() override;

         /**
          * @brief Activate output of parity splitting in power output
          */
         void showParity();

         /**
          * @brief Compute power for scalar field
          */
         void compute(Transform::TransformCoordinatorType& coord) override;

         /**
          * @brief Requires heavy calculation?
          */
         virtual bool isHeavy() const override; 
         
      protected:
         /**
          * @brief Write content
          */
         virtual void writeContent() override;

         /**
          * @brief Data ordering is m slowest
          */
         bool mHasMOrdering;

         /*
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
          * @brief Storage for the scalar power
          */
         Matrix mPower;

         /**
          * @brief Reset power storage
          */
         virtual void resetPower();

         /**
          * @brief Store power
          */
         virtual void storePower(const int n, const int l, const int m, const MHDFloat power);
   };

   inline bool ISphericalScalarRSpectrumWriter::isHeavy() const
   {
      return true;
   }

}
}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALSCALARRARSPECTRUMWRITER_HPP
