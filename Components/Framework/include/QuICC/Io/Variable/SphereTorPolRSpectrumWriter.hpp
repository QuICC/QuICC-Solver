/**
 * @file SphereTorPolRSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics radial power spectrum
 * calculation for a Toroidal/Poloidal field in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERETORPOLRSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERETORPOLRSPECTRUMWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Variable/ISphericalTorPolRSpectrumWriter.hpp"
#include "QuICC/Resolutions/Resolution.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

/**
 * @brief Implementation of the ASCII spherical harmonics radial power spectrum
 * calculation for a Toroidal/Poloidal field in a sphere
 */
class SphereTorPolRSpectrumWriter : public ISphericalTorPolRSpectrumWriter
{
public:
   /**
    * @brief Constructor
    *
    * @param prefix Prefix to use for file name
    * @param type Type of the file (typically scheme name)
    */
   SphereTorPolRSpectrumWriter(const std::string& prefix,
      const std::string& type);

   /**
    * @brief Destructor
    */
   virtual ~SphereTorPolRSpectrumWriter() = default;

   /**
    * @brief Initialise the operator, transform and file
    */
   virtual void init();

protected:
private:
};

/// Typedef for a shared pointer of a HDF5 state file writer
typedef std::shared_ptr<SphereTorPolRSpectrumWriter>
   SharedSphereTorPolRSpectrumWriter;

} // namespace Variable
} // namespace Io
} // namespace QuICC

#endif // QUICC_IO_VARIABLE_SPHERETORPOLRSPECTRUMWRITER_HPP
