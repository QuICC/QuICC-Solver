/**
 * @file ShellTorPolMSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics M energy spectrum
 * calculation for a Toroidal/Poloidal field in a spherical shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLTORPOLMSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLTORPOLMSPECTRUMWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Variable/ISphericalTorPolMSpectrumWriter.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

/**
 * @brief Implementation of the ASCII spherical harmonics M energy spectrum
 * calculation for a Toroidal/Poloidal field in a spherical shell
 */
class ShellTorPolMSpectrumWriter : public ISphericalTorPolMSpectrumWriter
{
public:
   /**
    * @brief Constructor
    *
    * @param prefix Prefix to use for file name
    * @param type Type of the file (typically scheme name)
    */
   ShellTorPolMSpectrumWriter(const std::string& prefix,
      const std::string& type,
      std::vector<std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction>> pF = {});

   /**
    * @brief Destructor
    */
   virtual ~ShellTorPolMSpectrumWriter() = default;

   /**
    * @brief Initialise the operator, transform and file
    */
   virtual void init();

protected:
private:
};

/// Typedef for a shared pointer of a HDF5 state file writer
typedef std::shared_ptr<ShellTorPolMSpectrumWriter>
   SharedShellTorPolMSpectrumWriter;

} // namespace Variable
} // namespace Io
} // namespace QuICC

#endif // QUICC_IO_VARIABLE_SHELLTORPOLMSPECTRUMWRITER_HPP
