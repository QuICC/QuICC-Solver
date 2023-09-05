/**
 * @file ISphericalTorPolModeSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics mode energy spectrum calculation for toroidal/poloidal field in a spherical geometry
 */

// System includes
//
#include <iomanip>
#include <stdexcept>

// Project includes
//
#include "QuICC/Io/Variable/ISphericalTorPolModeSpectrumWriter.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Io/Variable/Tags/Spectrum.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ISphericalTorPolModeSpectrumWriter::ISphericalTorPolModeSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolEnergyBaseWriter(prefix + Tags::Spectrum::MODEBASENAME, Tags::Spectrum::EXTENSION, prefix + Tags::Spectrum::HEADER, type, Tags::Spectrum::VERSION, Dimensions::Space::SPECTRAL, OVERWRITE), mTorEnergy(0,0), mPolEnergy(0,0)
   {
   }

   void ISphericalTorPolModeSpectrumWriter::init()
   {
      // Resize storage for spectra
      this->mTorEnergy = Matrix::Zero(this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL),this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL));
      this->mPolEnergy = Matrix::Zero(this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL),this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL));

      ISphericalTorPolEnergyBaseWriter::init();
   }

   void ISphericalTorPolModeSpectrumWriter::resetEnergy()
   {
      this->mTorEnergy.setZero();
      this->mPolEnergy.setZero();
   }

   void ISphericalTorPolModeSpectrumWriter::storeQEnergy(const int l, const int m, const MHDFloat energy)
   {
      if(m <= l)
      {
         this->mPolEnergy(l,m) += energy;
      }
   }

   void ISphericalTorPolModeSpectrumWriter::storeSEnergy(const int l, const int m, const MHDFloat energy)
   {
      if(m <= l)
      {
         this->mPolEnergy(l,m) += energy;
      }
   }

   void ISphericalTorPolModeSpectrumWriter::storeTEnergy(const int l, const int m, const MHDFloat energy)
   {
      if(m <= l)
      {
         this->mTorEnergy(l,m) += energy;
      }
   }

   void ISphericalTorPolModeSpectrumWriter::writeContent()
   {
      // Normalize by the volume
      this->mTorEnergy /= 2.0*this->mVolume;
      this->mPolEnergy /= 2.0*this->mVolume;

      // Create file
      this->preWrite();

      // Get the "global" energy spectrum from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, mTorEnergy.data(), mTorEnergy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, mPolEnergy.data(), mPolEnergy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      using Tools::Formatter::ioFW;
      using Tools::Formatter::ioIW;
      int ioPrec = 14;

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::scientific;
         this->mFile << "# time: " << std::setprecision(ioPrec);
         this->mFile << ioFW(ioPrec) << this->mTime << std::endl;

         this->mFile << "# energy: " << std::setprecision(ioPrec);
         MHDFloat tT = this->mTorEnergy.sum();
         MHDFloat tP = this->mPolEnergy.sum();
         MHDFloat t = tT + tP;
         this->mFile << ioFW(ioPrec) << t << "\t";
         this->mFile << ioFW(ioPrec) << tT << "\t";
         this->mFile << ioFW(ioPrec) << tP << std::endl;

         this->mFile << std::left << ioIW() << "l" << "\t" << "m" << "\t";
         this->mFile << ioFW(ioPrec) << "total" << "\t";
         this->mFile << ioFW(ioPrec) << "toroidal" << "\t";
         this->mFile << ioFW(ioPrec) << "poloidal" << std::endl;

         // Total
         assert(this->mTorEnergy.rows() == this->mPolEnergy.rows());
         for(int l = 0; l < this->mTorEnergy.rows(); l++)
         {
            for(int m = 0; m <= l; m++)
            {
               const MHDFloat eT = this->mTorEnergy(l,m);
               const MHDFloat eP = this->mPolEnergy(l,m);
               this->mFile << std::left << ioIW() << l << "\t" << m << "\t" << std::setprecision(ioPrec);
               this->mFile << ioFW(ioPrec) << eT + eP << "\t";
               this->mFile << ioFW(ioPrec) << eT << "\t";
               this->mFile << ioFW(ioPrec) << eP;
               this->mFile << std::endl;
            }
         }
      }

      // Close file
      this->postWrite();

      // Abort if energy spectrum contains NaN
      if(std::isnan(this->mTorEnergy.sum()) || std::isnan(this->mPolEnergy.sum()))
      {
         QuICCEnv().abort("Toroidal/Poloidal mode energy spectrum is NaN!");
      }
   }

}
}
}
