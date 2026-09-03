/**
 * @file SphereScalarNCoeffsPowerWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics L power spectrum
 * calculation for a scalar field in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERESCALARNCOEFFSPOWERWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERESCALARNCOEFFSPOWERWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Io/Variable/ISphericalScalarNCoeffsPowerWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

/**
 * @brief Implementation of the ASCII coefficients power spectrum
 * calculation for a scalar field in a sphere
 */
class SphereScalarNCoeffsPowerWriter : public ISphericalScalarNCoeffsPowerWriter
{
public:
   /**
    * @brief Constructor
    *
    * @param prefix Prefix to use for file name
    * @param type Type of the file (typically scheme name)
    */
   SphereScalarNCoeffsPowerWriter(const std::string& prefix,
      const std::string& type);

   /**
    * @brief Destructor
    */
   virtual ~SphereScalarNCoeffsPowerWriter() = default;

   /**
    * @brief Initialise the operator, transform and file
    */
   virtual void init();

protected:
private:
};

/// Typedef for a shared pointer of a HDF5 state file writer
typedef std::shared_ptr<SphereScalarNCoeffsPowerWriter>
   SharedSphereScalarNCoeffsPowerWriter;

} // namespace Variable
} // namespace Io
} // namespace QuICC

#endif // QUICC_IO_VARIABLE_SPHERESCALARNCOEFFSPOWERWRITER_HPP
