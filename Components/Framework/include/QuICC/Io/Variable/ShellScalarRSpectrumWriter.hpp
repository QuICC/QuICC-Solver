/**
 * @file ShellScalarRSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics radial power spectrum
 * calculation for a scalar field in a spherical shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLSCALARRSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLSCALARRSPECTRUMWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Io/Variable/ISphericalScalarRSpectrumWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

/**
 * @brief Implementation of the ASCII spherical harmonics radial power spectrum
 * calculation for a scalar field in a spherical shell
 */
class ShellScalarRSpectrumWriter : public ISphericalScalarRSpectrumWriter
{
public:
   /**
    * @brief Constructor
    *
    * @param prefix Prefix to use for file name
    * @param type Type of the file (typically scheme name)
    */
   ShellScalarRSpectrumWriter(const std::string& prefix,
      const std::string& type);

   /**
    * @brief Destructor
    */
   virtual ~ShellScalarRSpectrumWriter() = default;

   /**
    * @brief Initialise the operator, transform and file
    */
   virtual void init();

protected:
private:
};

/// Typedef for a shared pointer of a HDF5 state file writer
typedef std::shared_ptr<ShellScalarRSpectrumWriter>
   SharedShellScalarRSpectrumWriter;

} // namespace Variable
} // namespace Io
} // namespace QuICC

#endif // QUICC_IO_VARIABLE_SHELLSCALARRSPECTRUMWRITER_HPP
