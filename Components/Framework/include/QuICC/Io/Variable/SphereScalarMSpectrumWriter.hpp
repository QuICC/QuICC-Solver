/**
 * @file SphereScalarMSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics M energy spectrum
 * calculation for a scalar field in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERESCALARMSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERESCALARMSPECTRUMWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Io/Variable/ISphericalScalarMSpectrumWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

/**
 * @brief Implementation of the ASCII spherical harmonics M energy spectrum
 * calculation for a scalar field in a sphere
 */
class SphereScalarMSpectrumWriter : public ISphericalScalarMSpectrumWriter
{
public:
   /**
    * @brief Constructor
    *
    * @param prefix Prefix to use for file name
    * @param type Type of the file (typically scheme name)
    */
   SphereScalarMSpectrumWriter(const std::string& prefix,
      const std::string& type);

   /**
    * @brief Destructor
    */
   virtual ~SphereScalarMSpectrumWriter() = default;

   /**
    * @brief Initialise the operator, transform and file
    */
   virtual void init();

protected:
private:
};

/// Typedef for a shared pointer of a HDF5 state file writer
typedef std::shared_ptr<SphereScalarMSpectrumWriter>
   SharedSphereScalarMSpectrumWriter;

} // namespace Variable
} // namespace Io
} // namespace QuICC

#endif // QUICC_IO_VARIABLE_SPHERESCALARMSPECTRUMWRITER_HPP
