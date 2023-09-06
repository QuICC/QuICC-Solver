/** 
 * @file ICartesian1DScalarEnergyWriter.hpp
 * @brief Implementation of the ASCII Chebyshev energy calculation for a scalar field in a plane layer
 */

#ifndef QUICC_IO_VARIABLE_ICARTESIAN1DSCALARENERGYWRITER_HPP
#define QUICC_IO_VARIABLE_ICARTESIAN1DSCALARENERGYWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Io/Variable/ICartesian1DScalarEnergyBaseWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII Chebyshev energy calculation for a scalar field in a plane layer 
    */
   class ICartesian1DScalarEnergyWriter: public ICartesian1DScalarEnergyBaseWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ICartesian1DScalarEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ICartesian1DScalarEnergyWriter() = default;
         
      protected:
         /**
          * @brief Write content
          */
         virtual void writeContent();

      private:
         /**
          * @brief Storage for the scalar energy
          */
         Array mEnergy;

         /**
          * @brief Reset energy storage
          */
         virtual void resetEnergy();

         /**
          * @brief Store energy
          *
          * @param kx      Fourier mode in X
          * @param ky      Fourier mode in Y
          * @param energy  Energy of mode
          */
         virtual void storeEnergy(const int kx, const int ky, const MHDFloat energy);
   };

} // Variable
} // Io
} // QuICC

#endif // QUICC_IO_VARIABLE_ICARTESIAN1DSCALARENERGYWRITER_HPP
