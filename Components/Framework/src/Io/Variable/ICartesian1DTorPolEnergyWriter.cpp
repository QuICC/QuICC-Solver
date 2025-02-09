/**
 * @file ICartesian1DTorPolEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII Chebyshev energy calculation for toroidal/poloidal field in a plane layer
 */

// System includes
//
#include <iomanip>
#include <stdexcept>

// Project includes
//
#include "QuICC/Io/Variable/ICartesian1DTorPolEnergyWriter.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Io/Variable/Tags/Energy.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ICartesian1DTorPolEnergyWriter::ICartesian1DTorPolEnergyWriter(const std::string& prefix, const std::string& type)
      : ICartesian1DTorPolEnergyBaseWriter(prefix + Tags::Energy::BASENAME, Tags::Energy::EXTENSION, prefix + Tags::Energy::HEADER, type, Tags::Energy::VERSION, Dimensions::Space::SPECTRAL), mTorEnergy(2), mPolEnergy(2)
   {
      this->mTorEnergy.setConstant(-1);
      this->mPolEnergy.setConstant(-1);
   }

   void ICartesian1DTorPolEnergyWriter::resetEnergy()
   {
      this->mTorEnergy.setZero();
      this->mPolEnergy.setZero();
   }

   void ICartesian1DTorPolEnergyWriter::storeQEnergy(const int kx, const int ky, const MHDFloat energy)
   {
      if(kx == 0 && ky == 0)
      {
         this->mPolEnergy(0) += energy;
      } else
      {
         this->mPolEnergy(1) += energy;
      }
   }

   void ICartesian1DTorPolEnergyWriter::storeSEnergy(const int kx, const int ky, const MHDFloat energy)
   {
      if(kx == 0 && ky == 0)
      {
         this->mPolEnergy(0) += energy;
      } else
      {
         this->mPolEnergy(1) += energy;
      }
   }

   void ICartesian1DTorPolEnergyWriter::storeTEnergy(const int kx, const int ky, const MHDFloat energy)
   {
      if(kx == 0 && ky == 0)
      {
         this->mTorEnergy(0) += energy;
      } else
      {
         this->mTorEnergy(1) += energy;
      }
   }

   void ICartesian1DTorPolEnergyWriter::writeContent()
   {
      // Normalize by the volume
      this->mTorEnergy /= 2.0*this->mVolume;
      this->mPolEnergy /= 2.0*this->mVolume;

      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         Array energy(4);

         energy(0) = this->mTorEnergy(0);
         energy(1) = this->mTorEnergy(1);
         energy(2) = this->mPolEnergy(0);
         energy(3) = this->mPolEnergy(1);

         MPI_Allreduce(MPI_IN_PLACE, energy.data(), energy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mTorEnergy(0) = energy(0);
         this->mTorEnergy(1) = energy(1);
         this->mPolEnergy(0) = energy(2);
         this->mPolEnergy(1) = energy(3);
      #endif //QUICC_MPI

      using Tools::Formatter::ioFW;
      int ioPrec = 14;

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::scientific;
         // Total
         this->mFile << std::setprecision(ioPrec) << ioFW(ioPrec) << this->mTime << "\t" << ioFW(ioPrec) << this->mTorEnergy.sum() + this->mPolEnergy.sum();
         if(this->mShowParity)
         {
            this->mFile << std::setprecision(ioPrec) << "\t" << ioFW(ioPrec) << this->mTorEnergy(0) + this->mPolEnergy(0) << "\t" << ioFW(ioPrec) << this->mTorEnergy(1) + this->mPolEnergy(1);
         }

         // Toroidal
         this->mFile << std::setprecision(ioPrec) << "\t" << ioFW(ioPrec) << this->mTorEnergy.sum();
         if(this->mShowParity)
         {
            this->mFile << std::setprecision(ioPrec) << "\t" << ioFW(ioPrec) << this->mTorEnergy(0) << "\t" << ioFW(ioPrec) << this->mTorEnergy(1);
         }

         // Poloidal
         this->mFile << std::setprecision(ioPrec) << "\t" << ioFW(ioPrec) << this->mPolEnergy.sum();
         if(this->mShowParity)
         {
            this->mFile << std::setprecision(ioPrec) << "\t" << ioFW(ioPrec) << this->mPolEnergy(0) << "\t" << ioFW(ioPrec) << this->mPolEnergy(1);
         }

         // End line
         this->mFile << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mTorEnergy.sum()) || std::isnan(this->mPolEnergy.sum()))
      {
         QuICCEnv().abort("Toroidal/Poloidal energy is NaN!");
      }
   }

} // Variable
} // Io
} // QuICC
