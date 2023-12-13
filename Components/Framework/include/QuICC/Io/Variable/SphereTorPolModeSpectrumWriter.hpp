/**
 * @file SphereTorPolModeSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics mode energy spectrum
 * calculation for a Toroidal/Poloidal field in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERETORPOLMODESPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERETORPOLMODESPECTRUMWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Io/Variable/ISphericalTorPolModeSpectrumWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

/**
 * @brief Implementation of the ASCII spherical harmonics mode energy spectrum
 * calculation for a Toroidal/Poloidal field in a sphere
 */
class SphereTorPolModeSpectrumWriter : public ISphericalTorPolModeSpectrumWriter
{
public:
   /**
    * @brief Constructor
    *
    * @param prefix Prefix to use for file name
    * @param type Type of the file (typically scheme name)
    */
   SphereTorPolModeSpectrumWriter(const std::string& prefix,
      const std::string& type);

   /**
    * @brief Destructor
    */
   ~SphereTorPolModeSpectrumWriter() = default;

   /**
    * @brief Initialise the operator, transform and file
    */
   void init() final;

protected:
private:
};

/// Typedef for a shared pointer of a sphere Tor/Pol mode enery spectrum writer
typedef std::shared_ptr<SphereTorPolModeSpectrumWriter>
   SharedSphereTorPolModeSpectrumWriter;

} // namespace Variable
} // namespace Io
} // namespace QuICC

#endif // QUICC_IO_VARIABLE_SPHERETORPOLMODESPECTRUMWRITER_HPP
