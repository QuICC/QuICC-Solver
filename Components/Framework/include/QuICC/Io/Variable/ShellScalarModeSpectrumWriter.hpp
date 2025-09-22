/**
 * @file ShellScalarModeSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics mode energy spectrum
 * calculation for a scalar field in a spherical shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLSCALARMODESPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLSCALARMODESPECTRUMWRITER_HPP

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
 * calculation for a scalar field in a spherical shell
 */
class ShellScalarModeSpectrumWriter : public ISphericalScalarModeSpectrumWriter
{
public:
   /**
    * @brief Constructor
    *
    * @param prefix Prefix to use for file name
    * @param type Type of the file (typically scheme name)
    */
   ShellScalarModeSpectrumWriter(const std::string& prefix,
      const std::string& type);

   /**
    * @brief Destructor
    */
   virtual ~ShellScalarModeSpectrumWriter() = default;

   /**
    * @brief Initialise the operator, transform and file
    */
   virtual void init();

protected:
private:
};

/// Typedef for a shared pointer of a HDF5 state file writer
typedef std::shared_ptr<ShellScalarModeSpectrumWriter>
   SharedShellScalarModeSpectrumWriter;

} // namespace Variable
} // namespace Io
} // namespace QuICC

#endif // QUICC_IO_VARIABLE_SHELLSCALARMODESPECTRUMWRITER_HPP
