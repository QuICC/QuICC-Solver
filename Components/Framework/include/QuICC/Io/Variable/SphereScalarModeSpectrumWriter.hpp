/**
 * @file SphereScalarModeSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics mode energy spectrum
 * calculation for a scalar field in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERESCALARMODESPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERESCALARMODESPECTRUMWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Io/Variable/ISphericalScalarModeSpectrumWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

/**
 * @brief Implementation of the ASCII spherical harmonics mode energy spectrum
 * calculation for a scalar field in a sphere
 */
class SphereScalarModeSpectrumWriter : public ISphericalScalarModeSpectrumWriter
{
public:
   /**
    * @brief Constructor
    *
    * @param prefix Prefix to use for file name
    * @param type Type of the file (typically scheme name)
    */
   SphereScalarModeSpectrumWriter(const std::string& prefix,
      const std::string& type);

   /**
    * @brief Destructor
    */
   ~SphereScalarModeSpectrumWriter() = default;

   /**
    * @brief Initialise the operator, transform and file
    */
   void init() final;

protected:
private:
};

/// Typedef for a shared pointer of a sphere scalar spherical harmonic mode
/// energy spectrum writer
typedef std::shared_ptr<SphereScalarModeSpectrumWriter>
   SharedSphereScalarModeSpectrumWriter;

} // namespace Variable
} // namespace Io
} // namespace QuICC

#endif // QUICC_IO_VARIABLE_SPHERESCALARMODESPECTRUMWRITER_HPP
