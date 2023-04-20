/** 
 * @file ICartesian1DTorPolEnergyBaseWriter.hpp
 * @brief Implementation of the ASCII Chebyshev energy calculation for a Toroidal/Poloidal field in a plane layer
 */

#ifndef QUICC_IO_VARIABLE_ICARTESIAN1DTORPOLENERGYWBASERITER_HPP
#define QUICC_IO_VARIABLE_ICARTESIAN1DTORPOLENERGYWBASERITER_HPP

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
    * @brief Implementation of the ASCII Chebyshev energy calculation for a Toroidal/Poloidal field in a plane layer
    */
   class ICartesian1DTorPolEnergyBaseWriter: public IVariableAsciiWriter
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
         ICartesian1DTorPolEnergyBaseWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const IAsciiWriter::WriteMode mode = IAsciiWriter::EXTEND);

         /**
          * @brief Destructor
          */
         virtual ~ICartesian1DTorPolEnergyBaseWriter() = default;

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
          * @brief Store energy from Q component
          *
          * @param kx      Fourier mode in X
          * @param ky      Fourier mode in Y
          * @param energy  Energy of mode
          */
         virtual void storeQEnergy(const int kx, const int ky, const MHDFloat energy) = 0;

         /**
          * @brief Store energy from S component
          *
          * @param kx      Fourier mode in X
          * @param ky      Fourier mode in Y
          * @param energy  Energy of mode
          */
         virtual void storeSEnergy(const int kx, const int ky, const MHDFloat energy) = 0;

         /**
          * @brief Store energy from T component
          *
          * @param kx      Fourier mode in X
          * @param ky      Fourier mode in Y
          * @param energy  Energy of mode
          */
         virtual void storeTEnergy(const int kx, const int ky, const MHDFloat energy) = 0;
   };

   inline bool ICartesian1DTorPolEnergyBaseWriter::isHeavy() const
   {
      return true;
   }

} // Variable
} // Io
} // QuICC

#endif // QUICC_IO_VARIABLE_ICARTESIAN1DTORPOLENERGYBASEWRITER_HPP
