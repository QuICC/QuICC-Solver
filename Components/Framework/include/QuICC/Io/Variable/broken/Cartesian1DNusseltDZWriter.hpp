/** 
 * @file Cartesian1DNusseltDZWriter.hpp
 * @brief Implementation of the ASCII Nusselt number writer through the Z boundary extracted from temperature field
 */

#ifndef QUICC_IO_VARIABLE_CARTESIAN1DNUSSELTDZWRITER_HPP
#define QUICC_IO_VARIABLE_CARTESIAN1DNUSSELTDZWRITER_HPP

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
    * @brief Implementation of the ASCII Nusselt number writer through the Z boundary extracted from temperature field
    */
   class Cartesian1DNusseltDZWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DNusseltDZWriter(std::string type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DNusseltDZWriter();

         /**
          * @brief Initialise the operator and file
          */
         virtual void init();

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:

      private:
         /**
          * @brief Nusselt calculation operator
          */
         SparseMatrix   mNusseltOp;

   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<Cartesian1DNusseltDZWriter> SharedCartesian1DNusseltDZWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_CARTESIAN1DNUSSELTDZWRITER_HPP
