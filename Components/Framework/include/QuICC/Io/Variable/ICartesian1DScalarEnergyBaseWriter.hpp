/** 
 * @file ICartesian1DScalarEnergyBaseWriter.hpp
 * @brief Implementation of the ASCII Chebyshev energy calculation for a scalar field in a plane layer
 */

#ifndef QUICC_IO_VARIABLE_ICARTESIAN1DSCALARENERGYBASEWRITER_HPP
#define QUICC_IO_VARIABLE_ICARTESIAN1DSCALARENERGYBASEWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Io/Variable/IVariableAsciiWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII Chebyshev energy calculation for a scalar field in a plane layer
    */
   class ICartesian1DScalarEnergyBaseWriter: public IVariableAsciiWriter
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
         ICartesian1DScalarEnergyBaseWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const WriteMode mode = EXTEND);

         /**
          * @brief Destructor
          */
         virtual ~ICartesian1DScalarEnergyBaseWriter() = default;

         /**
          * @brief Activate output of parity splitting in energy output
          */
         void showParity();

         /**
          * @brief Compute energy for scalar field
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
         void prepareInput(Transform::TransformCoordinatorType& coord);

         /*
          * @brief Cartesian1D volume to normalize energy to energy density
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
          * @brief Store energy
          *
          * @param kx      Fourier mode in X
          * @param ky      Fourier mode in Y
          * @param energy  Energy of mode
          */
         virtual void storeEnergy(const int kx, const int ky, const MHDFloat energy) = 0;
   };

   inline bool ICartesian1DScalarEnergyBaseWriter::isHeavy() const
   {
      return true;
   }

} // Variable
} // Io
} // QuICC

#endif // QUICC_IO_VARIABLE_ICARTESIAN1DSCALARENERGYBASEWRITER_HPP
