/** 
 * @file ShellNusseltWriter.hpp
 * @brief Implementation of the Nusselt number in a spherical shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLNUSSELTWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLNUSSELTWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Io/Variable/IVariableAsciiWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the Nusselt number in a spherical shell
    */
   class ShellNusseltWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ShellNusseltWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ShellNusseltWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();

         /**
          * @brief Requires heavy calculation?
          */
         virtual bool isHeavy() const; 
         
      protected:
         /**
          * @brief Write State to file
          */
         virtual void writeContent();

         /**
          * @brief Data ordering is m slowest
          */
         bool mHasMOrdering;

      private:
         /**
          * @brief Nusselt number
          */
         Array mNusselt;

         /*
          * @brief Heat flux from background profile
          */
         Array mBackground;

         /**
          * @brief Origin projector
          */
         Matrix mBoundary;
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<ShellNusseltWriter> SharedShellNusseltWriter;

   inline bool ShellNusseltWriter::isHeavy() const
   {
      return false;
   }

}
}
}

#endif // QUICC_IO_VARIABLE_SHELLNUSSELTWRITER_HPP
