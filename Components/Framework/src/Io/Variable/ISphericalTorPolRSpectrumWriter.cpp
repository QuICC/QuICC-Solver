/**
 * @file ISphericalTorPolRSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics power spectrum calculation for toroidal/poloidal field in a spherical geometry
 */

// Configuration includes
//

// System includes
//
#include <iomanip>

// External includes
//

// Class include
//
#include "QuICC/Io/Variable/ISphericalTorPolRSpectrumWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/Transform/Reductor/RadialPower.hpp"
#include "QuICC/Transform/Reductor/RadialPowerOverr1.hpp"
#include "QuICC/Transform/Reductor/RadialPowerOverr1D1R1.hpp"
#include "QuICC/Io/Variable/Tags/Spectrum.hpp"

namespace QuICC {

namespace Io {

namespace Variable {
   ISphericalTorPolRSpectrumWriter::ISphericalTorPolRSpectrumWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + Tags::Spectrum::RBASENAME, Tags::Spectrum::EXTENSION, prefix + Tags::Spectrum::HEADER, type, Tags::Spectrum::VERSION, Dimensions::Space::SPECTRAL, OVERWRITE), mHasMOrdering(false), mVolume(std::numeric_limits<MHDFloat>::quiet_NaN()), mShowParity(false), mGrid(0), mTorPower(0,0), mPolPower(0,0)
   {
   }

   ISphericalTorPolRSpectrumWriter::~ISphericalTorPolRSpectrumWriter()
   {
   }

   void ISphericalTorPolRSpectrumWriter::init()
   {
      // Resize storage for spectra
      this->mTorPower = Matrix::Zero(this->res().sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL), this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL));
      this->mPolPower = Matrix::Zero(this->res().sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL), this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL));
   }

   void ISphericalTorPolRSpectrumWriter::showParity()
   {
      this->mShowParity = true;
   }

   void ISphericalTorPolRSpectrumWriter::prepareInput(const FieldComponents::Spectral::Id sId, Transform::TransformCoordinatorType& coord)
   {
      // Initialize grid
      if(this->mGrid.size() == 0)
      {
         this->mGrid = coord.transform1D().meshGrid();
      }

      // get iterator to field
      vector_iterator vIt;
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);
      auto&& field = vRange.first->second;
      assert(std::visit([&](auto&& p)->bool{return (p->dom(0).res().sim().ss().spectral().ONE() == FieldComponents::Spectral::TOR);}, field));
      assert(std::visit([&](auto&& p)->bool{return (p->dom(0).res().sim().ss().spectral().TWO() == FieldComponents::Spectral::POL);}, field));

      constexpr auto TId = Dimensions::Transform::TRA1D;
      const int packs = 1;
      coord.communicator().converter<TId>().setupCommunication(packs, TransformDirection::BACKWARD);

      coord.communicator().converter<TId>().prepareBackwardReceive();

      // Dealias variable data
      std::visit(
            [&](auto&& p)
            {
               coord.communicator().transferForward(Dimensions::Transform::SPECTRAL, p->rDom(0).rTotal().rComp(sId), false);
            },
            field);

      coord.communicator().converter<TId>().initiateForwardSend();
   }

   void ISphericalTorPolRSpectrumWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      constexpr auto TId = Dimensions::Transform::TRA1D;
      Matrix spectrum;

      // Prepare spectral data for transform
      this->prepareInput(FieldComponents::Spectral::TOR, coord);

      // Recover dealiased BWD data
      auto pInVarTor = coord.ss().bwdPtr(TId);
      coord.communicator().receiveBackward(TId, pInVarTor);

      const auto& tRes = *this->res().cpu()->dim(TId);

      // Size of spectrum
      auto size = std::visit([](auto&& p)->std::pair<int,int>{return std::make_pair(0,p->data().cols());}, pInVarTor);
      size.first = tRes.dim<Dimensions::Data::DATF1D>();

      // Compute power reduction
      spectrum.resize(size.first, size.second);
      std::visit(
            [&](auto&& p)
            {
               coord.transform1D().reduce(spectrum, p->data(), Transform::Reductor::RadialPower::id());
            },
            pInVarTor);

      this->resetPower();

      MHDFloat lfactor = 0.0;
      MHDFloat factor = 1.0;
      int idx = 0;
      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // m = 0, no factor of two
            int m_ = tRes.idx<Dimensions::Data::DAT3D>(k);
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
               lfactor = l_*(l_+1.0);

               for(int i = 0; i < tRes.dim<Dimensions::Data::DATF1D>(j, k); i++)
               {
                  this->storeTPower(i, l_, m_, factor*lfactor*spectrum(i, idx));
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
            lfactor = l_*(l_+1.0);
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

               for(int i = 0; i < tRes.dim<Dimensions::Data::DATF1D>(j, k); i++)
               {
                  this->storeTPower(i, l_, m_, factor*lfactor*spectrum(i, idx));
               }
               idx += 1;
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<TId>().freeBwd(pInVarTor);

      // Prepare spectral data for transform
      this->prepareInput(FieldComponents::Spectral::POL, coord);

      // Recover dealiased BWD data
      auto pInVarPolQ = coord.ss().bwdPtr(TId);
      coord.communicator().receiveBackward(TId, pInVarPolQ);

      // Compute power reduction
      spectrum.setZero();
      std::visit(
            [&](auto&& p)
            {
               coord.transform1D().reduce(spectrum, p->data(), Transform::Reductor::RadialPowerOverr1::id());
            },
            pInVarPolQ);

      // Compute power in Q component of QST decomposition
      idx = 0;
      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // m = 0, no factor of two
            int m_ = tRes.idx<Dimensions::Data::DAT3D>(k);
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
               lfactor = std::pow(l_*(l_+1.0),2);

               for(int i = 0; i < tRes.dim<Dimensions::Data::DATF1D>(j, k); i++)
               {
                  this->storeQPower(i, l_, m_, factor*lfactor*spectrum(i, idx));
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
            lfactor = std::pow(l_*(l_+1.0),2);
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

               for(int i = 0; i < tRes.dim<Dimensions::Data::DATF1D>(j, k); i++)
               {
                  this->storeQPower(i, l_, m_, factor*lfactor*spectrum(i, idx));
               }
               idx += 1;
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<TId>().freeBwd(pInVarPolQ);

      // Prepare spectral data for transform
      this->prepareInput(FieldComponents::Spectral::POL, coord);

      // Recover dealiased BWD data
      auto pInVarPolS = coord.ss().bwdPtr(TId);
      coord.communicator().receiveBackward(TId, pInVarPolS);

      // Compute power reduction
      spectrum.setZero();
      std::visit(
            [&](auto&& p)
            {
               coord.transform1D().reduce(spectrum, p->data(), Transform::Reductor::RadialPowerOverr1D1R1::id());
            },
            pInVarPolS);

      // Compute power in S component of QST decomposition
      idx = 0;
      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // m = 0, no factor of two
            int m_ = tRes.idx<Dimensions::Data::DAT3D>(k);
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
               lfactor = l_*(l_+1.0);

               for(int i = 0; i < tRes.dim<Dimensions::Data::DATF1D>(j, k); i++)
               {
                  this->storeSPower(i, l_, m_, factor*lfactor*spectrum(i, idx));
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
            lfactor = l_*(l_+1.0);
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

               for(int i = 0; i < tRes.dim<Dimensions::Data::DATF1D>(j, k); i++)
               {
                  this->storeSPower(i, l_, m_, factor*lfactor*spectrum(i, idx));
               }
               idx += 1;
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<TId>().freeBwd(pInVarPolS);
   }

   void ISphericalTorPolRSpectrumWriter::resetPower()
   {
      this->mTorPower.setZero();
      this->mPolPower.setZero();
   }

   void ISphericalTorPolRSpectrumWriter::storeQPower(const int n, const int l, const int m, const MHDFloat power)
   {
      this->mPolPower(n, l) += power;
   }

   void ISphericalTorPolRSpectrumWriter::storeSPower(const int n, const int l, const int m, const MHDFloat power)
   {
      this->mPolPower(n, l) += power;
   }

   void ISphericalTorPolRSpectrumWriter::storeTPower(const int n, const int l, const int m, const MHDFloat power)
   {
      this->mTorPower(n, l) += power;
   }

   void ISphericalTorPolRSpectrumWriter::writeContent()
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
         Array w = this->mGrid;
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
            this->mFile << std::left << ioIW() << this->mGrid(i) << "\t" << std::setprecision(ioPrec);
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
         QuICCEnv().abort("Toroidal/Poloidal N power spectrum is NaN!");
      }
   }

}
}
}
