/** 
 * @file ShellAngularMomentumWriter.hpp
 * @brief Implementation of the angular momentum in a shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLANGULARMOMENTUMWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLANGULARMOMENTUMWRITER_HPP

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
    * @brief Implementation of the angular momentum number in a shell
    */
   class ShellAngularMomentumWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ShellAngularMomentumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ShellAngularMomentumWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();

         /**
          * @brief Compute angular momentum
          */
         void compute(Transform::TransformCoordinatorType& coord);

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
          * @brief m = 0 mode is present
          */
         bool mHasM0;

         /**
          * @brief m = 1 mode is present
          */
         bool mHasM1;

         /**
          * @brief j index of m = 0 mode
          */
         int mM0j;

         /**
          * @brief j index of m = 0 mode
          */
         int mM0k;

         /**
          * @brief j index of m = 1 mode
          */
         int mM1j;

         /**
          * @brief j index of m = 1 mode
          */
         int mM1k;

         /**
          * @brief Angular momentum conservation operator
          */
         Matrix mOp;

         /**
          * @Brief X,Y,Z components of angular momentum
          */
         Array mMomentum;
   };

   /// Typedef for a shared pointer to and AngularMomentumWriter
   typedef std::shared_ptr<ShellAngularMomentumWriter> SharedShellAngularMomentumWriter;

   inline bool ShellAngularMomentumWriter::isHeavy() const
   {
      return true;
   }

}
}
}

#endif // QUICC_IO_VARIABLE_SHELLANGULARMOMENTUMWRITER_HPP
