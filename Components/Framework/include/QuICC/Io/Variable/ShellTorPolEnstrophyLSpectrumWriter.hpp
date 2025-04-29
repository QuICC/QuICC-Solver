/**
 * @file ShellTorPolEnstrophyLSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics enstrophy L Spectrum
 * calculation for a Toroidal/Poloidal field in a spherical shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLTORPOLENSTROPHYLSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLTORPOLENSTROPHYLSPECTRUMWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Io/Variable/ISphericalTorPolEnstrophyLSpectrumWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

/**
 * @brief Implementation of the ASCII spherical harmonics enstrophy L spectrum
 * calculation for a Toroidal/Poloidal field in a spherical shell
 */
class ShellTorPolEnstrophyLSpectrumWriter
    : public ISphericalTorPolEnstrophyLSpectrumWriter
{
public:
   /**
    * @brief Constructor
    *
    * @param prefix Prefix to use for file name
    * @param type Type of the file (typically scheme name)
    */
   ShellTorPolEnstrophyLSpectrumWriter(const std::string& prefix,
      const std::string& type);

   /**
    * @brief Destructor
    */
   virtual ~ShellTorPolEnstrophyLSpectrumWriter() = default;

   /**
    * @brief Initialise the operator, transform and file
    */
   virtual void init();

protected:
private:
};

/// Typedef for a shared pointer
typedef std::shared_ptr<ShellTorPolEnstrophyLSpectrumWriter>
   SharedShellTorPolEnstrophyLSpectrumWriter;

} // namespace Variable
} // namespace Io
} // namespace QuICC

#endif // QUICC_IO_VARIABLE_SHELLTORPOLENSTROPHYLSPECTRUMWRITER_HPP
