/**
 * @file ISphericalScalarEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for scalar field in a spherical geometry
 */

// System includes
//
#include <iomanip>
#include <stdexcept>

// Project includes
//
#include "QuICC/Io/Variable/ISphericalScalarEnergyWriter.hpp"
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/Tags/Energy.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ISphericalScalarEnergyWriter::ISphericalScalarEnergyWriter(const std::string& prefix, const std::string& type)
      : ISphericalScalarEnergyBaseWriter(prefix + Tags::Energy::BASENAME, Tags::Energy::EXTENSION, prefix + Tags::Energy::HEADER, type, Tags::Energy::VERSION, Dimensions::Space::SPECTRAL, ISphericalScalarEnergyBaseWriter::EXTEND), mEnergy(2)
   {
      this->mEnergy.setConstant(-1);
   }

   void ISphericalScalarEnergyWriter::resetEnergy()
   {
      this->mEnergy.setZero();
   }

   void ISphericalScalarEnergyWriter::storeEnergy(const int l, const int m, const MHDFloat energy)
   {
      if((l - m)%2 == 0)
      {
         this->mEnergy(0) += energy;
      } else
      {
         this->mEnergy(1) += energy;
      }
   }

   void ISphericalScalarEnergyWriter::writeContent()
   {
      // Normalize by the volume
      this->mEnergy /= this->mVolume;

      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mEnergy.data(), this->mEnergy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      using Tools::Formatter::ioFW;
      int ioPrec = 14;

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::scientific;
         this->mFile << std::setprecision(ioPrec) << ioFW(ioPrec) << this->mTime << "\t" << ioFW(ioPrec) << this->mEnergy.sum();
         if(this->mShowParity)
         {
            this->mFile << std::setprecision(ioPrec) << "\t" << ioFW(ioPrec) << this->mEnergy(0) << "\t" << ioFW(ioPrec) << this->mEnergy(1);
         }

         // End line
         this->mFile << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mEnergy.sum()))
      {
         QuICCEnv().abort("Scalar energy is NaN!");
      }
   }

} // Variable
} // Io
} // QuICC
