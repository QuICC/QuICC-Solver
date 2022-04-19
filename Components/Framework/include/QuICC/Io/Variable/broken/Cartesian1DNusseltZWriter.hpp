/** 
 * @file Cartesian1DNusseltZWriter.hpp
 * @brief Implementation of the ASCII Nusselt number writer through Z boundary
 */

#ifndef QUICC_IO_VARIABLE_CARTESIAN1DNUSSELTZWRITER_HPP
#define QUICC_IO_VARIABLE_CARTESIAN1DNUSSELTZWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Io/Variable/IVariableAsciiWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII Nusselt number writer
    */
   class Cartesian1DNusseltZWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DNusseltZWriter(std::string type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DNusseltZWriter();

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:

      private:

   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<Cartesian1DNusseltZWriter> SharedCartesian1DNusseltZWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_CARTESIAN1DNUSSELTZWRITER_HPP
