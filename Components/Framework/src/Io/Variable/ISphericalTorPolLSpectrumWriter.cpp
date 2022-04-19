/**
 * @file ISphericalTorPolLSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics L energy spectrum calculation for toroidal/poloidal field in a spherical geometry
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
#include "QuICC/Io/Variable/ISphericalTorPolLSpectrumWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Io/Variable/Tags/Spectrum.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ISphericalTorPolLSpectrumWriter::ISphericalTorPolLSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolEnergyBaseWriter(prefix + Tags::Spectrum::LBASENAME, Tags::Spectrum::EXTENSION, prefix + Tags::Spectrum::HEADER, type, Tags::Spectrum::VERSION, Dimensions::Space::SPECTRAL, OVERWRITE), mTorEnergy(0), mPolEnergy(0)
   {
   }

   ISphericalTorPolLSpectrumWriter::~ISphericalTorPolLSpectrumWriter()
   {
   }

   void ISphericalTorPolLSpectrumWriter::init()
   {
      // Resize storage for spectra
      this->mTorEnergy = Array::Zero(this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL));
      this->mPolEnergy = Array::Zero(this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL));

      ISphericalTorPolEnergyBaseWriter::init();
   }

   void ISphericalTorPolLSpectrumWriter::initializeEnergy()
   {
      this->mTorEnergy.setZero();
      this->mPolEnergy.setZero();
   }

   void ISphericalTorPolLSpectrumWriter::storeQEnergy(const int l, const int m, const MHDFloat energy)
   {
      this->mPolEnergy(l) += energy;
   }

   void ISphericalTorPolLSpectrumWriter::storeSEnergy(const int l, const int m, const MHDFloat energy)
   {
      this->mPolEnergy(l) += energy;
   }

   void ISphericalTorPolLSpectrumWriter::storeTEnergy(const int l, const int m, const MHDFloat energy)
   {
      this->mTorEnergy(l) += energy;
   }

   void ISphericalTorPolLSpectrumWriter::writeContent()
   {
      // Normalize by the volume
      this->mTorEnergy /= 2.0*this->mVolume;
      this->mPolEnergy /= 2.0*this->mVolume;

      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         Array energy(2*this->mTorEnergy.size());

         energy.segment(0,this->mTorEnergy.size()) = this->mTorEnergy;
         energy.segment(this->mTorEnergy.size(),this->mPolEnergy.size()) = this->mPolEnergy;

         MPI_Allreduce(MPI_IN_PLACE, energy.data(), energy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mTorEnergy = energy.segment(0,this->mTorEnergy.size());
         this->mPolEnergy = energy.segment(this->mTorEnergy.size(),this->mPolEnergy.size());
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

         this->mFile << "# SD: standard deviation" << std::endl;
         this->mFile << "# characteristic L: " << std::setprecision(ioPrec);
         Array w = Array::LinSpaced(this->mTorEnergy.size(), 0, this->mTorEnergy.size()-1);
         MHDFloat eT = (w.array()*this->mTorEnergy.array()).sum()/tT;
         MHDFloat eP = (w.array()*this->mPolEnergy.array()).sum()/tP;
         MHDFloat e = (w.array()*(this->mTorEnergy.array() + this->mPolEnergy.array())).sum()/t;
         MHDFloat vT = ((w.array() - eT).array().pow(2)*this->mTorEnergy.array()).sum()/tT;
         MHDFloat vP = ((w.array() - eP).array().pow(2)*this->mPolEnergy.array()).sum()/tP;
         MHDFloat v = ((w.array() - e).array().pow(2)*(this->mTorEnergy.array() + this->mPolEnergy.array())).sum()/t;
         this->mFile << std::setprecision(2) << std::fixed;
         this->mFile << e << " (SD=" << std::sqrt(v) << ")" << std::left << "\t";
         this->mFile << eT << " (SD=" << std::sqrt(vT) << ")" << std::left << "\t";
         this->mFile << eP << " (SD=" << std::sqrt(vP) << ")" << std::left;
         this->mFile << std::endl;
         this->mFile << std::scientific;

         this->mFile << std::left << ioIW() << "# l" << "\t";
         this->mFile << ioFW(ioPrec) << "total" << "\t";
         this->mFile << ioFW(ioPrec) << "toroidal" << "\t";
         this->mFile << ioFW(ioPrec) << "poloidal" << std::endl;

         // Total
         for(int i = 0; i < this->mTorEnergy.size(); i++)
         {
            this->mFile << std::left << ioIW() << i << "\t" << std::setprecision(ioPrec);
            this->mFile << ioFW(ioPrec) << this->mTorEnergy(i) + this->mPolEnergy(i) << "\t";
            this->mFile << ioFW(ioPrec) << this->mTorEnergy(i) << "\t";
            this->mFile << ioFW(ioPrec) << this->mPolEnergy(i);
            this->mFile << std::endl;
         }
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mTorEnergy.sum()) || std::isnan(this->mPolEnergy.sum()))
      {
         QuICCEnv().abort(99);

         throw std::logic_error("Toroidal/Poloidal L energy spectrum is NaN!");
      }
   }

}
}
}
