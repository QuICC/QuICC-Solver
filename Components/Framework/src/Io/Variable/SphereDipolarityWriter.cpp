/**
 * @file SphereDipolarityWriter.cpp
 * @brief Source of the implementation of the ASCII dipolarity in a sphere
 */

// System includes
//
#include <iomanip>
#include <stdexcept>

// Project includes
//
#include "QuICC/Io/Variable/SphereDipolarityWriter.hpp"
#include "Environment/QuICCEnv.hpp"
#include "Types/Math.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Io/Variable/Tags/Dipolarity.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Value.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   SphereDipolarityWriter::SphereDipolarityWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + Tags::Dipolarity::BASENAME, Tags::Dipolarity::EXTENSION, prefix + Tags::Dipolarity::HEADER, type, Tags::Dipolarity::VERSION, Dimensions::Space::SPECTRAL, EXTEND), mHasMOrdering(false), mAxialDipole(0.0), mNonAxialDipole(0.0)
   {
   }

   void SphereDipolarityWriter::init()
   {
      this->mHasMOrdering = this->res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123);
      const auto& tRes = *this->res().cpu()->dim(Dimensions::Transform::SPECTRAL);

      // Compute boundary operators
      int nN = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nL = this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);
      this->mCmbSpectrum.resize(nL);

      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int m_ = tRes.idx<Dimensions::Data::DAT3D>(k);
            for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);
               if(this->mValue.count(l_) == 0)
               {
                  this->mValue.try_emplace(l_, Array());
               }
            }
         }
      } else
      {
         // Loop over harmonic degree l
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l_ = tRes.idx<Dimensions::Data::DAT3D>(k);
            this->mValue.try_emplace(l_, Array());
         }
      }

      Polynomial::Worland::WorlandBase wb;
      for (auto& [l, op] : this->mValue)
      {
         auto a = wb.alpha(l);
         auto db = wb.dBeta();      
         SparseSM::Worland::Boundary::Value bc(a, db, l);
         op = bc.compute(nN-1).cast<MHDFloat>();
      }

      IVariableAsciiWriter::init();
   }

   void SphereDipolarityWriter::prepareInput(const FieldComponents::Spectral::Id sId, Transform::TransformCoordinatorType& coord)
   {
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

   void SphereDipolarityWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      constexpr auto TId = Dimensions::Transform::TRA1D;
      MatrixZ spectrum;
      this->resetEnergy();

      const auto& tRes = *this->res().cpu()->dim(TId);

      // Process ununused toroidal component
      this->prepareInput(FieldComponents::Spectral::TOR, coord);
      auto pInVarTor = coord.ss().bwdPtr(TId);
      coord.communicator().receiveBackward(TId, pInVarTor);
      coord.communicator().storage<TId>().freeBwd(pInVarTor);


      // Prepare spectral data for transform
      this->prepareInput(FieldComponents::Spectral::POL, coord);

      // Recover dealiased BWD data
      auto pInVarPolQ = coord.ss().bwdPtr(TId);
      coord.communicator().receiveBackward(TId, pInVarPolQ);

      // Compute energy reduction
      spectrum.resize(std::visit([](auto&& p)->int{return p->data().cols();}, pInVarPolQ), 1);

      // Compute CMB magnetic field spectrum B^2
      MHDFloat factor, lfactor;
      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // m = 0, no factor of two
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
               lfactor = l_ * l_ * (2.0 * l_ + 1.0) * (l_ + 1.0);

               std::visit(
                     [&](auto&& p)
                     {
                        spectrum = this->mValue.at(l_).transpose() * p->profile(j,k); 
                     },
                     pInVarPolQ);

               this->mCmbSpectrum(l_) += factor * lfactor * (spectrum.array().abs2()).real()(0);

               if (l_ == 1 and m_ == 0)
               {
                  this->mAxialDipole = (l_+1.)*l_*spectrum.real()(0);
               }
               if (l_ == 1 and m_ == 1)
               {
                  this->mNonAxialDipole = (l_+1.)*l_*spectrum(0);
               }
            }
         }
      } else
      {
         // Loop over harmonic degree l
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l_ = tRes.idx<Dimensions::Data::DAT3D>(k);
            lfactor = l_ * l_ * (2.0 * l_ + 1.0) * (l_ + 1.0);
            const auto& op = this->mValue.at(l_);
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

               std::visit(
                     [&](auto&& p)
                     {
                        spectrum = this->mValue.at(l_).transpose() * p->profile(j,k); 
                     },
                     pInVarPolQ);

               this->mCmbSpectrum(l_) += factor * lfactor * (spectrum.array().abs2()).real()(0);

               if (l_ == 1 and m_ == 0)
               {
                  this->mAxialDipole = (l_+1.)*l_*spectrum.real()(0);
               }
               if (l_ == 1 and m_ == 1)
               {
                  this->mNonAxialDipole = (l_+1.)*l_*spectrum(0);
               }
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<TId>().freeBwd(pInVarPolQ);

      // Process unused poloidal field
      this->prepareInput(FieldComponents::Spectral::POL, coord);
      auto pInVarPolS = coord.ss().bwdPtr(TId);
      coord.communicator().receiveBackward(TId, pInVarPolS);
      coord.communicator().storage<TId>().freeBwd(pInVarPolS);
   }

   void SphereDipolarityWriter::writeContent()
   {
      // Create file
      this->preWrite();

      // Get the "global" value from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mCmbSpectrum.data(), this->mCmbSpectrum.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &this->mAxialDipole, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &this->mNonAxialDipole, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      using Tools::Formatter::ioFW;
      int ioPrec = 14;

      // Compute dipolarity
      this->mDipolarity = std::sqrt(this->mCmbSpectrum(1) / this->mCmbSpectrum.topRows(13).sum());

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::scientific;
         this->mFile << std::setprecision(ioPrec) << ioFW(ioPrec) << this->mTime << "\t" << ioFW(ioPrec) << this->mDipolarity << "\t" << this->mAxialDipole << "\t" << this->mNonAxialDipole.real() << "\t" << this->mNonAxialDipole.imag();
         this->mFile << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if is NaN
      if(std::isnan(this->mCmbSpectrum.sum()))
      {
         QuICCEnv().abort("Sphere dipolarity momentum is NaN!");
      }
   }

   void SphereDipolarityWriter::resetEnergy()
   {
      this->mCmbSpectrum.setZero();
      this->mAxialDipole = 0.0;
      this->mNonAxialDipole = 0.0;
   }

} // namespace Io
} // namespace Variable
} // namespace QuICC
