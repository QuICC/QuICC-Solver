/** 
 * @file Cartesian1DScalarEnergyWriter.hpp
 * @brief Implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a scalar field
 */

#ifndef QUICC_IO_VARIABLE_CARTESIAN1DSCALARENERGYWRITER_HPP
#define QUICC_IO_VARIABLE_CARTESIAN1DSCALARENERGYWRITER_HPP

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
    * @brief Implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a scalar field
    */
   class Cartesian1DScalarEnergyWriter: public IVariableAsciiWriter
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
         virtual ~Cartesian1DScalarEnergyWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();

         /**
          * @brief Compute energy for scalar field
          */
         void compute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief Write State to file
          */
         virtual void write();

         /**
          * @brief Requires heavy calculation?
          */
         virtual bool isHeavy() const; 
         
      protected:

      private:
         /*
          * @brief Cartesian box volume to normalize energy to energy density
          */
         MHDFloat mVolume;

         /**
          * @brief Storage for the scalar energy
          */
         Array mEnergy;

   };

   inline bool Cartesian1DScalarEnergyWriter::isHeavy() const
   {
      return true;
   }

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<Cartesian1DScalarEnergyWriter> SharedCartesian1DScalarEnergyWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_CARTESIAN1DSCALARENERGYWRITER_HPP
