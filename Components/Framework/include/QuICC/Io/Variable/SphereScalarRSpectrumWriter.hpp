/** 
 * @file SphereScalarRSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics radial power spectrum calculation for a scalar field in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERESCALARRSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERESCALARRSPECTRUMWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Io/Variable/ISphericalScalarRSpectrumWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics radial power spectrum calculation for a scalar field in a sphere
    */
   class SphereScalarRSpectrumWriter: public ISphericalScalarRSpectrumWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         SphereScalarRSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~SphereScalarRSpectrumWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<SphereScalarRSpectrumWriter> SharedSphereScalarRSpectrumWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_SPHERESCALARRSPECTRUMWRITER_HPP
