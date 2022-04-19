/** 
 * @file SphereScalarNSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics L power spectrum calculation for a scalar field in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERESCALARNSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERESCALARNSPECTRUMWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Io/Variable/ISphericalScalarNSpectrumWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics L power spectrum calculation for a scalar field in a sphere
    */
   class SphereScalarNSpectrumWriter: public ISphericalScalarNSpectrumWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         SphereScalarNSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~SphereScalarNSpectrumWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<SphereScalarNSpectrumWriter> SharedSphereScalarNSpectrumWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_SPHERESCALARNSPECTRUMWRITER_HPP
