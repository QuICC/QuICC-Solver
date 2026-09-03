/**
 * @file SphereTorPolLSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics L energy spectrum
 * calculation for a Toroidal/Poloidal field in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERETORPOLLSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERETORPOLLSPECTRUMWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Variable/ISphericalTorPolLSpectrumWriter.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "DenseSM/Worland/RadialTorPolFunction.hpp"


namespace QuICC {

namespace Io {

namespace Variable {

/**
 * @brief Implementation of the ASCII spherical harmonics L energy spectrum
 * calculation for a Toroidal/Poloidal field in a sphere
 */
class SphereTorPolLSpectrumWriter : public ISphericalTorPolLSpectrumWriter
{
public:
   /**
    * @brief Constructor
    *
    * @param prefix Prefix to use for file name
    * @param type Type of the file (typically scheme name)
    */
   SphereTorPolLSpectrumWriter(const std::string& prefix,
      const std::string& type, std::vector<std::shared_ptr<QuICC::DenseSM::Worland::RadialTorPolFunction>> pF = {});

   /**
    * @brief Destructor
    */
   virtual ~SphereTorPolLSpectrumWriter() = default;

   /**
    * @brief Initialise the operator, transform and file
    */
   virtual void init();

protected:
private:
};

/// Typedef for a shared pointer of a HDF5 state file writer
typedef std::shared_ptr<SphereTorPolLSpectrumWriter>
   SharedSphereTorPolLSpectrumWriter;

} // namespace Variable
} // namespace Io
} // namespace QuICC

#endif // QUICC_IO_VARIABLE_SPHERETORPOLLSPECTRUMWRITER_HPP
