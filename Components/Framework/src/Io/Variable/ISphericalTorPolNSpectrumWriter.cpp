/**
 * @file ISphericalTorPolNSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics N power spectrum calculation for toroidal/poloidal field in a spherical geometry
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
#include "QuICC/Io/Variable/ISphericalTorPolNSpectrumWriter.hpp"

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

   ISphericalTorPolNSpectrumWriter::ISphericalTorPolNSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolPowerBaseWriter(prefix + Tags::Spectrum::NBASENAME, Tags::Spectrum::EXTENSION, prefix + Tags::Spectrum::HEADER, type, Tags::Spectrum::VERSION, Dimensions::Space::SPECTRAL, OVERWRITE), mTorPower(0,0), mPolPower(0,0)
   {
   }

   ISphericalTorPolNSpectrumWriter::~ISphericalTorPolNSpectrumWriter()
   {
   }

   void ISphericalTorPolNSpectrumWriter::init()
   {
      // Resize storage for spectra
      this->mTorPower = Matrix::Zero(this->res().sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL), this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL));
      this->mPolPower = Matrix::Zero(this->res().sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL), this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL));

      ISphericalTorPolPowerBaseWriter::init();
   }

   void ISphericalTorPolNSpectrumWriter::initializePower()
   {
      this->mTorPower.setZero();
      this->mPolPower.setZero();
   }

   void ISphericalTorPolNSpectrumWriter::storeQPower(const int n, const int l, const int m, const MHDFloat power)
   {
      this->mPolPower(n, l) += power;
   }

   void ISphericalTorPolNSpectrumWriter::storeSPower(const int n, const int l, const int m, const MHDFloat power)
   {
      this->mPolPower(n, l) += power;
   }

   void ISphericalTorPolNSpectrumWriter::storeTPower(const int n, const int l, const int m, const MHDFloat power)
   {
      this->mTorPower(n, l) += power;
   }

   void ISphericalTorPolNSpectrumWriter::writeContent()
   {
      // Normalize by the volume
      this->mTorPower /= 2.0*this->mVolume;
      this->mPolPower /= 2.0*this->mVolume;

      // Create file
      this->preWrite();

      // Get the "global" Kinetic power from MPI code
      #ifdef QUICC_MPI
         Matrix power(this->mTorPower.rows(),2*this->mTorPower.cols());

         power.leftCols(this->mTorPower.cols()) = this->mTorPower;
         power.rightCols(this->mTorPower.cols()) = this->mPolPower;

         MPI_Allreduce(MPI_IN_PLACE, power.data(), power.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mTorPower = power.leftCols(this->mTorPower.cols());
         this->mPolPower = power.rightCols(this->mTorPower.cols());
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
         MHDFloat tT = this->mTorPower.sum();
         MHDFloat tP = this->mPolPower.sum();
         MHDFloat t = tT + tP;
         this->mFile << ioFW(ioPrec) << t << "\t";
         this->mFile << ioFW(ioPrec) << tT << "\t";
         this->mFile << ioFW(ioPrec) << tP << std::endl;

         this->mFile << "# SD: standard deviation" << std::endl;
         this->mFile << "# PC: percent of total" << std::endl;
         this->mFile << "# characteristic N:" << std::endl;
         Array w = Array::LinSpaced(this->mTorPower.rows(), 0, this->mTorPower.rows()-1);
         for(int l = 1; l < this->mTorPower.cols(); l++)
         {
            MHDFloat pT = this->mTorPower.col(l).sum();
            MHDFloat pP = this->mPolPower.col(l).sum();
            MHDFloat p = pT + pP;
            MHDFloat eT = (w.array()*this->mTorPower.col(l).array()).sum()/pT;
            MHDFloat eP = (w.array()*this->mPolPower.col(l).array()).sum()/pP;
            MHDFloat e = (w.array()*(this->mTorPower.col(l).array() + this->mPolPower.col(l).array())).sum()/p;
            MHDFloat vT = ((w.array() - eT).array().pow(2)*this->mTorPower.col(l).array()).sum()/pT;
            MHDFloat vP = ((w.array() - eP).array().pow(2)*this->mPolPower.col(l).array()).sum()/pP;
            MHDFloat v = ((w.array() - e).array().pow(2)*(this->mTorPower.col(l).array() + this->mPolPower.col(l).array())).sum()/p;
            this->mFile << "# l = " << ioIW() << l << "\t" << std::setprecision(2) << std::fixed;
            this->mFile << e << " (SD=" << std::sqrt(v) << ",pc=" << std::setprecision(1) << std::right << std::setw(4) << (p/t)*100 << ")" << std::left << "\t";
            this->mFile << eT << " (SD=" << std::sqrt(vT) << ",pc=" << std::setprecision(1) << std::right << std::setw(4) << (pT/tT)*100 << ")" << std::left << "\t";
            this->mFile << eP << " (SD=" << std::sqrt(vP) << ",pc=" << std::setprecision(1) << std::right << std::setw(4) << (pP/tP)*100 << ")" << std::left;
            this->mFile << std::endl;
            this->mFile << std::scientific;
         }

         this->mFile << std::left << ioIW() << "# n" << "\t";
         this->mFile << ioFW(ioPrec) << "Total";
         this->mFile << std::endl;

         // Total spectrum
         for(int i = 0; i < this->mTorPower.rows(); i++)
         {
            this->mFile << std::left << ioIW() << i << "\t" << std::setprecision(ioPrec);
            for(int j = 0; j < this->mTorPower.cols(); j++)
            {
               this->mFile << ioFW(ioPrec) << this->mTorPower(i,j) + this->mPolPower(i,j) << "\t";
            }
            this->mFile << std::endl;
         }

         this->mFile << std::endl << std::endl;
         this->mFile << std::left << ioIW() << "# n" << "\t";
         this->mFile << ioFW(ioPrec) << "Toroidal";
         this->mFile << std::endl;

         // Toroidal spectrum
         for(int i = 0; i < this->mTorPower.rows(); i++)
         {
            this->mFile << std::left << ioIW() << i << "\t" << std::setprecision(ioPrec);
            for(int j = 0; j < this->mTorPower.cols(); j++)
            {
               this->mFile << ioFW(ioPrec) << this->mTorPower(i,j) << "\t";
            }
            this->mFile << std::endl;
         }

         this->mFile << std::endl << std::endl;
         this->mFile << std::left << ioIW() << "# n" << "\t";
         this->mFile << ioFW(ioPrec) << "Poloidal";
         this->mFile << std::endl;

         // Poloidal spectrum
         for(int i = 0; i < this->mPolPower.rows(); i++)
         {
            this->mFile << std::left << ioIW() << i << "\t" << std::setprecision(ioPrec);
            for(int j = 0; j < this->mPolPower.cols(); j++)
            {
               this->mFile << ioFW(ioPrec) << this->mPolPower(i,j) << "\t";
            }
            this->mFile << std::endl;
         }
      }

      // Close file
      this->postWrite();

      // Abort if kinetic power is NaN
      if(std::isnan(this->mTorPower.sum()) || std::isnan(this->mPolPower.sum()))
      {
         QuICCEnv().abort(99);

         throw std::logic_error("Toroidal/Poloidal N power spectrum is NaN!");
      }
   }

}
}
}
