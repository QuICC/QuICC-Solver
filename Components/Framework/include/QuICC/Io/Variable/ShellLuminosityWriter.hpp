/** 
 * @file ShellLuminosityWriter.hpp
 * @brief Implementation of the Luminosity in a spherical shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLLUMINOSITYWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLLUMINOSITYWRITER_HPP

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
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the Luminosity in a spherical shell
    */
   class ShellLuminosityWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ShellLuminosityWriter(const std::string& prefix, 
                               const std::string& type, 
                               std::vector<std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction>> pF); 

         /**
          * @brief Destructor
          */
         virtual ~ShellLuminosityWriter();

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

         /**
          * @brief Vector of shared pointers to radial profiles (e.g. density)
          */
         std::vector<std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction>> mpF;

         /**
          * @brief shared pointers to density*Temperature*kappa profile 
          */
         std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> mpRhoTempKappa;

         /**
          * @brief shared pointers to D1ConductiveEntropy profile 
          */
         std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> mpD1Sc;

      private:
         /**
          * @brief Luminosity
          */
         Array mLuminosity;

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
   typedef std::shared_ptr<ShellLuminosityWriter> SharedShellLuminosityWriter;

   inline bool ShellLuminosityWriter::isHeavy() const
   {
      return false;
   }

}
}
}

#endif // QUICC_IO_VARIABLE_SHELLLUMINOSITYWRITER_HPP
