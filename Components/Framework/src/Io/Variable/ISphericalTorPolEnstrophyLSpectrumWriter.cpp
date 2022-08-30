/**
 * @file ISphericalTorPolEnstrophyLSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics enstrophy L spectrum calculation for toroidal/poloidal field in a spherical geometry
 */

// Configuration includes
//

// System includes
//
#include <iomanip>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Io/Variable/ISphericalTorPolEnstrophyLSpectrumWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Io/Variable/Tags/EnstrophySpectrum.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ISphericalTorPolEnstrophyLSpectrumWriter::ISphericalTorPolEnstrophyLSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolEnstrophyBaseWriter(prefix + Tags::EnstrophySpectrum::BASENAME + Tags::EnstrophySpectrum::LBASENAME, Tags::EnstrophySpectrum::EXTENSION, prefix + Tags::EnstrophySpectrum::HEADER, type, Tags::EnstrophySpectrum::VERSION, Dimensions::Space::SPECTRAL), mTorEnstrophy(0), mPolEnstrophy(0)
   {
   }

   ISphericalTorPolEnstrophyLSpectrumWriter::~ISphericalTorPolEnstrophyLSpectrumWriter()
   {
   }

   void ISphericalTorPolEnstrophyLSpectrumWriter::init()
   {
      // Resize storage for spectra
      this->mTorEnstrophy = Array::Zero(this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL));
      this->mPolEnstrophy = Array::Zero(this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL));

      ISphericalTorPolEnstrophyBaseWriter::init();
   }

   void ISphericalTorPolEnstrophyLSpectrumWriter::initializeEnstrophy()
   {
      this->mTorEnstrophy.setZero();
      this->mPolEnstrophy.setZero();
   }

   void ISphericalTorPolEnstrophyLSpectrumWriter::storeQEnstrophy(const int l, const int m, const MHDFloat enstrophy)
   {
      this->mTorEnstrophy(l) += enstrophy;
   }

   void ISphericalTorPolEnstrophyLSpectrumWriter::storeSEnstrophy(const int l, const int m, const MHDFloat enstrophy)
   {
      this->mTorEnstrophy(l) += enstrophy;
   }

   void ISphericalTorPolEnstrophyLSpectrumWriter::storeTEnstrophy(const int l, const int m, const MHDFloat enstrophy)
   {
      this->mPolEnstrophy(l) += enstrophy;
   }

   void ISphericalTorPolEnstrophyLSpectrumWriter::writeContent()
   {
      // Normalize by the volume
      this->mTorEnstrophy /= this->mVolume;
      this->mPolEnstrophy /= this->mVolume;

      // Create file
      this->preWrite();

      // Get the "global" Kinetic enstrophy from MPI code
      #ifdef QUICC_MPI
         Array enstrophy(2*this->mTorEnstrophy.size());

         enstrophy.segment(0,this->mTorEnstrophy.size()) = this->mTorEnstrophy;
         enstrophy.segment(this->mTorEnstrophy.size(),this->mPolEnstrophy.size()) = this->mPolEnstrophy;

         MPI_Allreduce(MPI_IN_PLACE, enstrophy.data(), enstrophy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mTorEnstrophy = enstrophy.segment(0,this->mTorEnstrophy.size());
         this->mPolEnstrophy = enstrophy.segment(this->mTorEnstrophy.size(),this->mPolEnstrophy.size());
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

         this->mFile << "# enstrophy: " << std::setprecision(ioPrec);
         this->mFile << ioFW(ioPrec) << this->mTorEnstrophy.sum() + this->mPolEnstrophy.sum() << "\t";
         this->mFile << ioFW(ioPrec) << this->mTorEnstrophy.sum() << "\t";
         this->mFile << ioFW(ioPrec) << this->mPolEnstrophy.sum() << std::endl;

         this->mFile << std::left << ioIW() << "# l" << "\t";
         this->mFile << ioFW(ioPrec) << "total" << "\t";
         this->mFile << ioFW(ioPrec) << "toroidal" << "\t";
         this->mFile << ioFW(ioPrec) << "poloidal" << std::endl;

         // Total
         for(int i = 0; i < this->mTorEnstrophy.size(); i++)
         {
            this->mFile << std::left << ioIW() << i << "\t" << std::setprecision(ioPrec);
            this->mFile << ioFW(ioPrec) << this->mTorEnstrophy(i) + this->mPolEnstrophy(i) << "\t";
            this->mFile << ioFW(ioPrec) << this->mTorEnstrophy(i) << "\t";
            this->mFile << ioFW(ioPrec) << this->mPolEnstrophy(i) << std::endl;
         }
      }

      // Close file
      this->postWrite();

      // Abort if kinetic enstrophy is NaN
      if(std::isnan(this->mTorEnstrophy.sum()) || std::isnan(this->mPolEnstrophy.sum()))
      {
         QuICCEnv().abort("Toroidal/Poloidal enstrophy L spectrum is NaN!");
      }
   }

}
}
}
