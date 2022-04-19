/** 
 * @file Cartesian3DNusseltZWriter.hpp
 * @brief Implementation of the ASCII Nusselt number writer for a 3D box through Z
 */

#ifndef QUICC_IO_VARIABLE_CARTESIAN3DNUSSELTZWRITER_HPP
#define QUICC_IO_VARIABLE_CARTESIAN3DNUSSELTZWRITER_HPP

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
    * @brief Implementation of the ASCII Nusselt number writer for a 3D box through Z
    */
   class Cartesian3DNusseltZWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type Type of the file (typically scheme name)
          */
         Cartesian3DNusseltZWriter(std::string type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian3DNusseltZWriter();

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
   typedef std::shared_ptr<Cartesian3DNusseltZWriter> SharedCartesian3DNusseltZWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_CARTESIAN3DNUSSELTZWRITER_HPP
