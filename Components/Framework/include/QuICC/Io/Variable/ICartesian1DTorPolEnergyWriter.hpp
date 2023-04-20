/** 
 * @file ICartesian1DTorPolEnergyWriter.hpp
 * @brief Implementation of the ASCII Chebyshev energy calculation for a Toroidal/Poloidal field in a plane layer
 */

#ifndef QUICC_IO_VARIABLE_ICARTESIAN1DTORPOLENERGYWRITER_HPP
#define QUICC_IO_VARIABLE_ICARTESIAN1DTORPOLENERGYWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Io/Variable/ICartesian1DTorPolEnergyBaseWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII Chebyshev energy calculation for a Toroidal/Poloidal field in a plane layer
    */
   class ICartesian1DTorPolEnergyWriter: public ICartesian1DTorPolEnergyBaseWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ICartesian1DTorPolEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ICartesian1DTorPolEnergyWriter() = default;
         
      protected:
         /**
          * @brief Write content
          */
         virtual void writeContent();

      private:
         /**
          * @brief Storage for the Toroidal energy
          */
         Array mTorEnergy;

         /**
          * @brief Storage for the Poloidal energy
          */
         Array mPolEnergy;

         /**
          * @brief Reset energy storage
          */
         virtual void resetEnergy();

         /**
          * @brief Store energy from Q component
          *
          * @param kx      Fourier mode in X
          * @param ky      Fourier mode in Y
          * @param energy  Energy of mode
          */
         virtual void storeQEnergy(const int kx, const int ky, const MHDFloat energy);

         /**
          * @brief Store energy from S component
          *
          * @param kx      Fourier mode in X
          * @param ky      Fourier mode in Y
          * @param energy  Energy of mode
          */
         virtual void storeSEnergy(const int kx, const int ky, const MHDFloat energy);

         /**
          * @brief Store energy from T component
          *
          * @param kx      Fourier mode in X
          * @param ky      Fourier mode in Y
          * @param energy  Energy of mode
          */
         virtual void storeTEnergy(const int kx, const int ky, const MHDFloat energy);
   };

} // Variable
} // Io
} // QuICC

#endif // QUICC_IO_VARIABLE_ICARTESIAN1DTORPOLENERGYWRITER_HPP
