/**
 * @file ShellTorPolDissipationSpectraWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics dissipation calculation for a Toroidal/Poloidal field in a spherical shell
 * @author Nicolò Lardelli \<nicolo.lardelli@erdw.ethz.ch\>
 */
#ifndef QUICC_IO_VARIABLE_SHELLTORPOLDISSIPATIONSPECTRAWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLTORPOLDISSIPATIONSPECTRAWRITER_HPP

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
      class ShellTorPolDissipationSpectraWriter : public IVariableAsciiWriter {
      public:
          /**
           * @brief Constructor
           *
           * @param prefix Prefix to use for file name
           * @param type Type of the file (typically scheme name)
           */
          ShellTorPolDissipationSpectraWriter(const std::string &prefix, const std::string &type);

          /**
           * @brief Destructor
           */
          virtual ~ShellTorPolDissipationSpectraWriter();

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
          Matrix mTorDiss;

          /**
           * @brief Storage for the Poloidal energy
           */
          Matrix mPolDiss;

          /*
           * @brief Storage for the radial spectral decomposition, unused
           */
          Array mTorRadial;
          Array mPolRadial;
          /**
             * @brief Chebyshev operator to integrate in radius
             */
          SparseMatrix mIntgOp;

          /**
           * @brief Chebyshev operator for spherical integral in radius (include r^2 factor)
           */
          SparseMatrix mSphIntgOp;
      };

      inline bool ShellTorPolDissipationSpectraWriter::isHeavy() const {
         return true;
      }

      /// Typedef for a shared pointer of a HDF5 state file writer
      typedef std::shared_ptr<ShellTorPolDissipationSpectraWriter> SharedShellTorPolDissipationSpectraWriter;

   }
}
}

#endif //QUICC_IO_VARIABLE_SHELLTORPOLDISSIPATIONSPECTRAWRITER_HPP
