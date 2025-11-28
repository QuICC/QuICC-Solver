/**
 * @file ISphericalScalarNSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics L power spectrum
 * calculation for a scalar field in a spherical geometry
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALSCALARNSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALSCALARNSPECTRUMWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Variable/ISphericalScalarPowerBaseWriter.hpp"
#include "QuICC/Resolutions/Resolution.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

/**
 * @brief Implementation of the ASCII spherical harmonics L power spectrum
 * calculation for a scalar field in a spherical geometry
 */
class ISphericalScalarNSpectrumWriter : public ISphericalScalarPowerBaseWriter
{
public:
   /**
    * @brief Constructor
    *
    * @param prefix Prefix to use for file name
    * @param type Type of the file (typically scheme name)
    */
   ISphericalScalarNSpectrumWriter(const std::string& prefix,
      const std::string& type);

   /**
    * @brief Destructor
    */
   virtual ~ISphericalScalarNSpectrumWriter() = default;

   /**
    * @brief Initialise the operator, transform and file
    */
   virtual void init() override;

protected:
   /**
    * @brief Write content
    */
   virtual void writeContent() override;

private:
   /**
    * @brief Storage for the scalar power
    */
   Matrix mPower;

   /**
    * @brief Reset power storage
    */
   virtual void resetPower() override;

   /**
    * @brief Store power
    */
   virtual void storePower(const int n, const int l, const int m,
      const MHDFloat power) override;
};

} // namespace Variable
} // namespace Io
} // namespace QuICC

#endif // QUICC_IO_VARIABLE_ISPHERICALSCALARNSPECTRUMWRITER_HPP
