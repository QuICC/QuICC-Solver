/** 
 * @file Cartesian1DScalarEnergyWriter.hpp
 * @brief Implementation of the ASCII plane layer energy calculation for a scalar field
 */

#ifndef QUICC_IO_VARIABLE_CARTESIAN1DSCALARENERGYWRITER_HPP
#define QUICC_IO_VARIABLE_CARTESIAN1DSCALARENERGYWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Io/Variable/ICartesian1DScalarEnergyWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII plane layer energy calculation for a scalar field
    */
   class Cartesian1DScalarEnergyWriter: public ICartesian1DScalarEnergyWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DScalarEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DScalarEnergyWriter() = default;

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<Cartesian1DScalarEnergyWriter> SharedSCartesian1DScalarEnergyWriter;

} // Variable
} // Io
} // QuICC

#endif // QUICC_IO_VARIABLE_CARTESIAN1DSCALARENERGYWRITER_HPP
