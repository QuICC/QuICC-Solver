/**
 * @file ISphericalScalarModeSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics mode energy spectrum
 * calculation for a scalar field in a spherical geometry
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALSCALARMODESPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALSCALARMODESPECTRUMWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Io/Variable/ISphericalScalarEnergyBaseWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

/**
 * @brief Implementation of the ASCII spherical harmonics mode energy spectrum
 * calculation for a scalar field in a spherical geometry
 */
class ISphericalScalarModeSpectrumWriter
    : public ISphericalScalarEnergyBaseWriter
{
public:
   /**
    * @brief Constructor
    *
    * @param prefix Prefix to use for file name
    * @param type Type of the file (typically scheme name)
    */
   ISphericalScalarModeSpectrumWriter(const std::string& prefix,
      const std::string& type);

   /**
    * @brief Destructor
    */
   virtual ~ISphericalScalarModeSpectrumWriter() = default;

   /**
    * @brief Initialise the operator, transform and file
    */
   virtual void init();

protected:
   /**
    * @brief Write content
    */
   virtual void writeContent();

private:
   /**
    * @brief Storage for the scalar energy
    */
   Matrix mEnergy;

   /**
    * @brief Reset energy storage
    */
   virtual void resetEnergy();

   /**
    * @brief Store energy
    */
   virtual void storeEnergy(const int l, const int m, const MHDFloat energy);
};

} // namespace Variable
} // namespace Io
} // namespace QuICC

#endif // QUICC_IO_VARIABLE_ISPHERICALSCALARMODESPECTRUMWRITER_HPP
