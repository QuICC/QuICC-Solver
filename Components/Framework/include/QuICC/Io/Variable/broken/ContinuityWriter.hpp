/** 
 * @file ContinuityWriter.hpp
 * @brief Implementation of the ASCII maximal continuity writer
 */

#ifndef QUICC_IO_VARIABLE_CONTINUITYWRITER_HPP
#define QUICC_IO_VARIABLE_CONTINUITYWRITER_HPP

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
    * @brief Implementation of the ASCII maximal continuity writer
    */
   class ContinuityWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type Type of the file (typically scheme name)
          */
         ContinuityWriter(std::string type);

         /**
          * @brief Destructor
          */
         virtual ~ContinuityWriter();

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:

      private:

   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<ContinuityWriter> SharedContinuityWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_CONTINUITYWRITER_HPP
