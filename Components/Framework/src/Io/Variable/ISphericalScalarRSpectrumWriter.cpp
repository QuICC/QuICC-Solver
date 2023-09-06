/**
 * @file ISphericalScalarRSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics power on radial grid calculation for scalar field in a spherical geometry
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
#include "QuICC/Io/Variable/ISphericalScalarRSpectrumWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/Transform/Reductor/RadialPower.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/Tags/Spectrum.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ISphericalScalarRSpectrumWriter::ISphericalScalarRSpectrumWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + Tags::Spectrum::RBASENAME, Tags::Spectrum::EXTENSION, prefix + Tags::Spectrum::HEADER, type, Tags::Spectrum::VERSION, Dimensions::Space::SPECTRAL, OVERWRITE), mHasMOrdering(false), mVolume(std::numeric_limits<MHDFloat>::quiet_NaN()), mShowParity(false), mGrid(0), mPower(0,0)
   {
   }

   ISphericalScalarRSpectrumWriter::~ISphericalScalarRSpectrumWriter()
   {
   }

   void ISphericalScalarRSpectrumWriter::init()
   {
      this->mPower = Matrix::Zero(this->res().sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL), this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL));
   }

   void ISphericalScalarRSpectrumWriter::showParity()
   {
      this->mShowParity = true;
   }

   void ISphericalScalarRSpectrumWriter::resetPower()
   {
      this->mPower.setZero();
   }

   void ISphericalScalarRSpectrumWriter::storePower(const int n, const int l, const int m, const MHDFloat power)
   {
      this->mPower(n, l) += power;
   }

   void ISphericalScalarRSpectrumWriter::prepareInput(Transform::TransformCoordinatorType& coord)
   {
      // Initialize grid
      if(this->mGrid.size() == 0)
      {
         this->mGrid = coord.transform1D().meshGrid();
      }

      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);
      auto&& field = sRange.first->second;

      constexpr auto TId = Dimensions::Transform::TRA1D;
      const int packs = 1;
      coord.communicator().converter<TId>().setupCommunication(packs, TransformDirection::BACKWARD);

      coord.communicator().converter<TId>().prepareBackwardReceive();

      // Dealias variable data
      std::visit(
            [&](auto&& p)
            {
               coord.communicator().transferForward(Dimensions::Transform::SPECTRAL, p->rDom(0).rTotal(), false);
            },
            field);

      coord.communicator().converter<TId>().initiateForwardSend();
   }

   void ISphericalScalarRSpectrumWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      constexpr auto TId = Dimensions::Transform::TRA1D;

      // Prepare spectral data for transform
      this->prepareInput(coord);

      // Recover dealiased BWD data
      auto pInVar = coord.ss().bwdPtr(TId);
      coord.communicator().receiveBackward(TId, pInVar);

      const auto& tRes = *this->res().cpu()->dim(TId);

      // Size of spectrum
      auto size = std::visit([](auto&& p)->std::pair<int,int>{return std::make_pair(0,p->data().cols());}, pInVar);
      size.first = tRes.dim<Dimensions::Data::DATF1D>();

      // Compute power reduction
      Matrix spectrum(size.first, size.second);
      std::visit(
            [&](auto&& p)
            {
               coord.transform1D().reduce(spectrum, p->data(), Transform::Reductor::RadialPower::id());
            },
            pInVar);

      this->resetPower();

      MHDFloat factor = 1.0;
      int idx = 0;
      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int m_ = tRes.idx<Dimensions::Data::DAT3D>(k);
            // m = 0, no factor of two
            if(m_ == 0)
            {
               factor = 1.0;
            } else
            {
               factor = 2.0;
            }

            for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);

               for(int i = 0; i < tRes.dim<Dimensions::Data::DATF1D>(k); i++)
               {
                  this->storePower(i, l_, m_, factor*spectrum(i, idx));
               }
               idx += 1;
            }
         }
      } else
      {
         // Loop over harmonic degree l
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l_ = tRes.idx<Dimensions::Data::DAT3D>(k);

            // m = 0, no factor of two
            for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int m_ = tRes.idx<Dimensions::Data::DAT2D>(j,k);
               // m = 0, no factor of two
               if(m_ == 0)
               {
                  factor = 1.0;
               } else
               {
                  factor = 2.0;
               }

               for(int i = 0; i < tRes.dim<Dimensions::Data::DATF1D>(k); i++)
               {
                  this->storePower(i, l_, m_, factor*spectrum(i, idx));
               }
               idx += 1;
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<TId>().freeBwd(pInVar);
   }

   void ISphericalScalarRSpectrumWriter::writeContent()
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
         Array w = this->mGrid;
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
            this->mFile << std::left << ioIW() << this->mGrid(i) << "\t" << std::setprecision(ioPrec);
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
         QuICCEnv().abort("Scalar power is NaN!");
      }
   }

}
}
}
