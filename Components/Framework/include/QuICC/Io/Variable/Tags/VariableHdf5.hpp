/** 
 * @file VariableHdf5.hpp
 * @brief Definitions and names used by the spectral space data files
 */

#ifndef QUICC_IO_VARIABLE_TAGS_TAGS_VARIABLEHDF5_HPP
#define QUICC_IO_VARIABLE_TAGS_TAGS_VARIABLEHDF5_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Io {

namespace Variable {

namespace Tags {

   /**
    * @brief Implementation of HDF5File for a spectral data HDF5 file
    */
   class VariableHdf5
   {
      public:
         /**
          * @brief Truncation group tag name
          */
         static const std::string   TRUNCATION;

         /**
          * @brief Physical truncation group tag name
          */
         static const std::string   TRUNCPHYSICAL;

         /**
          * @brief Transform truncation group tag name
          */
         static const std::string   TRUNCTRANSFORM;

         /**
          * @brief Spectral truncation group tag name
          */
         static const std::string   TRUNCSPECTRAL;

         /**
          * @brief Truncation dimension tag name
          */
         static const std::string   TRUNCDIM;

         /**
          * @brief Physical parameters part for State file
          */
         static const std::string   PHYSICAL;

         /**
          * @brief Run parameters part for State file
          */
         static const std::string   RUN;

         /**
          * @brief Time tag for State file
          */
         static const std::string   RUNTIME;

         /**
          * @brief Timestep tag for State file
          */
         static const std::string   RUNSTEP;

      private:
         /**
         * @brief Empty destructor
         */
         VariableHdf5();

         /**
         * @brief Destructor
         */
         ~VariableHdf5();
   };
}
}
}
}

#endif // QUICC_IO_VARIABLE_TAGS_TAGS_VARIABLEHDF5_HPP
