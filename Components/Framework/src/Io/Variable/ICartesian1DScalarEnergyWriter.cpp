/**
 * @file ICartesian1DScalarEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII Chebyshev energy calculation for scalar field in a plane layer 
 */

// System includes
//
#include <iomanip>
#include <stdexcept>

// Project includes
//
#include "QuICC/Io/Variable/ICartesian1DScalarEnergyWriter.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/Tags/Energy.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ICartesian1DScalarEnergyWriter::ICartesian1DScalarEnergyWriter(const std::string& prefix, const std::string& type)
      : ICartesian1DScalarEnergyBaseWriter(prefix + Tags::Energy::BASENAME, Tags::Energy::EXTENSION, prefix + Tags::Energy::HEADER, type, Tags::Energy::VERSION, Dimensions::Space::SPECTRAL, ICartesian1DScalarEnergyBaseWriter::EXTEND), mEnergy(2)
   {
      this->mEnergy.setConstant(-1);
   }

   void ICartesian1DScalarEnergyWriter::resetEnergy()
   {
      this->mEnergy.setZero();
   }

   void ICartesian1DScalarEnergyWriter::storeEnergy(const int kx, const int ky, const MHDFloat energy)
   {
      if(kx == 0 && ky == 0)
      {
         this->mEnergy(0) += energy;
      } else
      {
         this->mEnergy(1) += energy;
      }
   }

   void ICartesian1DScalarEnergyWriter::writeContent()
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

}
}
}
