/** 
 * @file Cartesian1DMagneticEnergyWriter.hpp
 * @brief Implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a vector field (streamfunction formulation)
 */

#ifndef QUICC_IO_VARIABLE_CARTESIAN1DMAGNETICENERGYWRITER_HPP
#define QUICC_IO_VARIABLE_CARTESIAN1DMAGNETICENERGYWRITER_HPP

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
    * @brief Implementation of the ASCII Cartesian 1D (double periodic)energy calculation for a vector field (streamfunction formulation)
    */
   class Cartesian1DMagneticEnergyWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Enum type for different energies
          */
         enum EnergyTypeId {
            TOTAL,
            ZONAL_X,
            ZONAL_Y,
         };

         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DMagneticEnergyWriter(const std::string& prefix, const std::string& type, const bool hasZonalX = false, const bool hasZonalY = false);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DMagneticEnergyWriter();

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
         /**
          * @brief Compute energy for scalar field
          */
         void computeEnergy(Transform::TransformCoordinatorType& coord, const EnergyTypeId flag);

         /**
          * @brief Flag to compute zonal energy along X
          */
         bool mHasZonalX;

         /**
          * @brief Flag to compute zonal energy along Y
          */
         bool mHasZonalY;

         /**
          * @brief Cartesian box volume to normalize energy to energy density
          */
         MHDFloat mVolume;

         /**
          * @brief Storage for the X component energy
          */
         MHDFloat mXEnergy;

         /**
          * @brief Storage for the Y component energy
          */
         MHDFloat mYEnergy;

         /**
          * @brief Storage for the Z component energy
          */
         MHDFloat mZEnergy;

         /**
          * @brief Storage for the X component X zonal energy
          */
         MHDFloat mXZonalXEnergy;

         /**
          * @brief Storage for the Z component X zonal energy
          */
         MHDFloat mZZonalXEnergy;

         /**
          * @brief Storage for the Y component Y zonal energy
          */
         MHDFloat mYZonalYEnergy;

         /**
          * @brief Storage for the Z component Y zonal energy
          */
         MHDFloat mZZonalYEnergy;

         /**
          * @brief Chebyshev operator to integrate
          */
         SparseMatrix mIntgOp;
   };

   inline bool Cartesian1DMagneticEnergyWriter::isHeavy() const
   {
      return true;
   }

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<Cartesian1DMagneticEnergyWriter> SharedCartesian1DMagneticEnergyWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_CARTESIAN1DMAGNETICENERGYWRITER_HPP
