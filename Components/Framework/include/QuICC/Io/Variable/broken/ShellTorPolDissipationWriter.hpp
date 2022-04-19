/** 
 * @file ShellTorPolEnergyWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a spherical shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLTORPOLDISSIPATIONWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLTORPOLDISSIPATIONWRITER_HPP

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
       * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a spherical shell
       */
      class ShellTorPolDissipationWriter : public IVariableAsciiWriter {
      public:
          /**
           * @brief Constructor
           *
           * @param prefix Prefix to use for file name
           * @param type Type of the file (typically scheme name)
           */
          ShellTorPolDissipationWriter(const std::string &prefix, const std::string &type);

          /**
           * @brief Destructor
           */
          virtual ~ShellTorPolDissipationWriter();

          /**
           * @brief Initialise the operator, transform and file
           */
          virtual void init();

          /**
           * @brief Compute energy for scalar field
           */
          void compute(Transform::TransformCoordinatorType &coord);

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
           * @brief Spherical volume to normalize energy to energy density
           */
          MHDFloat mVolume;

          /**
           * @brief Storage for the Toroidal energy
           */
          MHDFloat mTorDiss;

          /**
           * @brief Storage for the Poloidal energy
           */
          MHDFloat mPolDiss;

          /**
           * @brief Storage for the centro-symmentric energy
           */
          MHDFloat mCentroSymDiss;

          /**
           * @brief Storage for the centro-antisymmetric energy
           */
          MHDFloat mCentroAntysymDiss;

          /**
          * @brief Storage for the equatorial-symmentric energy
          */
          MHDFloat mEquaSymDiss;

          /**
          * @brief Storage for the equatorial-antisymmetric energy
          */
          MHDFloat mEquaAntysymDiss;

          /**
             * @brief Chebyshev operator to integrate in radius
             */
          SparseMatrix mIntgOp;

          /**
           * @brief Chebyshev operator for spherical integral in radius (include r^2 factor)
           */
          SparseMatrix mSphIntgOp;
      };

      inline bool ShellTorPolDissipationWriter::isHeavy() const {
         return true;
      }

      /// Typedef for a shared pointer of a HDF5 state file writer
      typedef std::shared_ptr<ShellTorPolDissipationWriter> SharedShellTorPolDissipationWriter;

   }
}
}

#endif // QUICC_IO_VARIABLE_SHELLTORPOLDISSIPATIONWRITER_HPP
