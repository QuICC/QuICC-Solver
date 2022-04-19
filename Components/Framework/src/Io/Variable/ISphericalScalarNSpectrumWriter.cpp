/**
 * @file ISphericalScalarNSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics L power spectrum calculation for scalar field in a spherical geometry
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
#include "QuICC/Io/Variable/ISphericalScalarNSpectrumWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/Tags/Spectrum.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ISphericalScalarNSpectrumWriter::ISphericalScalarNSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalScalarPowerBaseWriter(prefix + Tags::Spectrum::NBASENAME, Tags::Spectrum::EXTENSION, prefix + Tags::Spectrum::HEADER, type, Tags::Spectrum::VERSION, Dimensions::Space::SPECTRAL, OVERWRITE), mPower(0,0)
   {
   }

   ISphericalScalarNSpectrumWriter::~ISphericalScalarNSpectrumWriter()
   {
   }

   void ISphericalScalarNSpectrumWriter::init()
   {
      this->mPower = Matrix::Zero(this->res().sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL), this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL));

      ISphericalScalarPowerBaseWriter::init();
   }

   void ISphericalScalarNSpectrumWriter::initializePower()
   {
      this->mPower.setZero();
   }

   void ISphericalScalarNSpectrumWriter::storePower(const int n, const int l, const int m, const MHDFloat power)
   {
      this->mPower(n, l) += power;
   }

   void ISphericalScalarNSpectrumWriter::writeContent()
   {
      // Normalize by the volume
      this->mPower /= this->mVolume;

      // Create file
      this->preWrite();

      // Get the "global" Kinetic power from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mPower.data(), this->mPower.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
         MHDFloat t = this->mPower.sum();
         this->mFile << ioFW(ioPrec) << t << std::endl;

         this->mFile << "# SD: standard deviation" << std::endl;
         this->mFile << "# PC: percent of total" << std::endl;
         this->mFile << "# characteristic N:" << std::endl;
         Array w = Array::LinSpaced(this->mPower.rows(), 0, this->mPower.rows()-1);
         for(int l = 0; l < this->mPower.cols(); l++)
         {
            MHDFloat p = this->mPower.col(l).sum();
            MHDFloat e = (w.array()*this->mPower.col(l).array()).sum()/p;
            MHDFloat v = ((w.array() - e).array().pow(2)*this->mPower.col(l).array()).sum()/p;
            this->mFile << "# l = " << ioIW() << l << "\t" << std::setprecision(2) << std::fixed;
            this->mFile << e << " (SD=" << std::sqrt(v) << ",pc=" << std::setprecision(1) << std::right << std::setw(4) << (p/t)*100 << ")" << std::left;
            this->mFile << std::endl;
            this->mFile << std::scientific;
         }

         this->mFile << std::left << ioIW() << "# n" << "\t";
         this->mFile << ioFW(ioPrec) << "Total";
         this->mFile << std::endl;

         // Spectrum
         for(int i = 0; i < this->mPower.rows(); i++)
         {
            this->mFile << std::left << ioIW() << i << "\t" << std::setprecision(ioPrec);
            for(int j = 0; j < this->mPower.cols(); j++)
            {
               this->mFile << ioFW(ioPrec) << this->mPower(i,j) << "\t";
            }
            this->mFile << std::endl;
         }

         // End line
         this->mFile << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic power is NaN
      if(std::isnan(this->mPower.sum()))
      {
         QuICCEnv().abort(99);

         throw std::logic_error("Scalar power is NaN!");
      }
   }

}
}
}
