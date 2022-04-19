/** 
 * @file SphereScalarLSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics L energy spectrum calculation for a scalar field in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERESCALARLSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERESCALARLSPECTRUMWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Io/Variable/ISphericalScalarLSpectrumWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII spherical harmonics L energy spectrum calculation for a scalar field in a sphere
    */
   class SphereScalarLSpectrumWriter: public ISphericalScalarLSpectrumWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         SphereScalarLSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~SphereScalarLSpectrumWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<SphereScalarLSpectrumWriter> SharedSphereScalarLSpectrumWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_SPHERESCALARLSPECTRUMWRITER_HPP
