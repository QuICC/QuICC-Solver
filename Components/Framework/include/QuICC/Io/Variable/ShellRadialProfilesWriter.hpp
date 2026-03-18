/** 
 * @file ShellRadialProfilesWriter.hpp
 * @brief Implementation of the Radial Profiles in a spherical shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLRADIALPROFILESWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLRADIALPROFILESWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Io/Variable/IVariableAsciiWriter.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"
#include "QuICC/Polynomial/Quadrature/ChebyshevRule.hpp"


namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the Radial Profiles in a spherical shell
    */
   class ShellRadialProfilesWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ShellRadialProfilesWriter(const std::string& prefix, 
                               const std::string& type, 
                               std::vector<std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction>> pF);

         /**
          * @brief Destructor
          */
         virtual ~ShellRadialProfilesWriter();

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

      private:

         /**
          * @brief Nusselt number
          */
         Internal::Array mGrid;

         /**
          * @brief Nusselt number
          */
         Internal::Matrix mProfiles;

         /**
          * @brief shared pointers to radial profiles vector
          */
         std::vector<std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction>> mPF;
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<ShellRadialProfilesWriter> SharedShellRadialProfilesWriter;

   inline bool ShellRadialProfilesWriter::isHeavy() const
   {
      return false;
   }

}
}
}

#endif // QUICC_IO_VARIABLE_SHELLRADIALPROFILESWRITER_HPP
