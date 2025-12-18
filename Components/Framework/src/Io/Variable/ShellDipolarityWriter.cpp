/**
 * @file ShellDipolarityWriter.cpp
 * @brief Source of the implementation of the ASCII dipolarity in a shell
 */

// System includes
//
#include <iomanip>
#include <stdexcept>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Variable/ShellDipolarityWriter.hpp"
#include "QuICC/Io/Variable/Tags/Dipolarity.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/ICondition.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/Value.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

ShellDipolarityWriter::ShellDipolarityWriter(const std::string& prefix,
   const std::string& type) :
    IVariableAsciiWriter(prefix + Tags::Dipolarity::BASENAME,
       Tags::Dipolarity::EXTENSION, prefix + Tags::Dipolarity::HEADER, type,
       Tags::Dipolarity::VERSION, Dimensions::Space::SPECTRAL, EXTEND),
    mHasMOrdering(false),
    mCmbNl(-1),
    mAxialDipole(0.0),
    mNonAxialDipole(0.0)
{}

void ShellDipolarityWriter::init()
{
   this->mHasMOrdering = this->res().sim().ss().has(
      SpatialScheme::Feature::TransformSpectralOrdering123);
   const auto& tRes = *this->res().cpu()->dim(Dimensions::Transform::TRA1D);

   // Compute boundary operators
   int nN = this->res().sim().dim(Dimensions::Simulation::SIM1D,
      Dimensions::Space::SPECTRAL);
   int nL = this->res().sim().dim(Dimensions::Simulation::SIM2D,
      Dimensions::Space::SPECTRAL);
   this->mCmbSpectrum.resize(nL);

   if (this->mHasMOrdering)
   {
      // Loop over harmonic order m
      for (int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
      {
         //int m_ = tRes.idx<Dimensions::Data::DAT3D>(k);
         for (int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
         {
            int l_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);
            if (this->mValue.count(l_) == 0)
            {
               this->mValue.try_emplace(l_, Array());
            }
         }
      }
   }
   else
   {
      // Loop over harmonic degree l
      for (int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
      {
         int l_ = tRes.idx<Dimensions::Data::DAT3D>(k);
         this->mValue.try_emplace(l_, Array());
      }
   }

   typedef SparseSM::Chebyshev::LinearMap::Boundary::ICondition::Position
      Position;
   MHDFloat ri =
      this->mPhysical.find(NonDimensional::Lower1d::id())->second->value();
   MHDFloat ro =
      this->mPhysical.find(NonDimensional::Upper1d::id())->second->value();

   for (auto& [l, op]: this->mValue)
   {
      SparseSM::Chebyshev::LinearMap::Boundary::Value bc(ri, ro, Position::TOP);
      op = bc.compute(nN - 1).cast<MHDFloat>();
   }

   if(this->mCmbNl < 0)
   {
      this->mCmbNl = this->mCmbSpectrum.size();
   }

   IVariableAsciiWriter::init();
}

void ShellDipolarityWriter::prepareInput(
   const FieldComponents::Spectral::Id sId,
   Transform::TransformCoordinatorType& coord)
{
   // get iterator to field
   vector_iterator vIt;
   vector_iterator_range vRange = this->vectorRange();
   assert(std::distance(vRange.first, vRange.second) == 1);
   auto&& field = vRange.first->second;
   assert(std::visit(
      [&](auto&& p) -> bool
      {
         return (p->dom(0).res().sim().ss().spectral().ONE() ==
                 FieldComponents::Spectral::TOR);
      },
      field));
   assert(std::visit(
      [&](auto&& p) -> bool
      {
         return (p->dom(0).res().sim().ss().spectral().TWO() ==
                 FieldComponents::Spectral::POL);
      },
      field));

   constexpr auto TId = Dimensions::Transform::TRA1D;
   const int packs = 1;
   coord.communicator().converter<TId>().setupCommunication(packs,
      TransformDirection::BACKWARD);

   coord.communicator().converter<TId>().prepareBackwardReceive();

   // Dealias variable data
   std::visit(
      [&](auto&& p)
      {
         coord.communicator().transferForward(Dimensions::Transform::SPECTRAL,
            p->rDom(0).rTotal().rComp(sId), false);
      },
      field);

   coord.communicator().converter<TId>().initiateForwardSend();
}

void ShellDipolarityWriter::compute(Transform::TransformCoordinatorType& coord)
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
   spectrum.resize(
      std::visit([](auto&& p) -> int { return p->data().cols(); }, pInVarPolQ),
      1);

   // Radial truncation
   int nN = this->res().sim().dim(Dimensions::Simulation::SIM1D,
      Dimensions::Space::SPECTRAL);

   // Compute CMB magnetic field spectrum B^2
   MHDFloat factor, lfactor;

   if (this->mHasMOrdering)
   {
      // Loop over harmonic order m
      for (int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
      {
         // m = 0, no factor of two
         int m_ = tRes.idx<Dimensions::Data::DAT3D>(k);
         // m = 0, no factor of two
         if (m_ == 0)
         {
            factor = 1.0;
         }
         else
         {
            factor = 2.0;
         }

         for (int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
         {
            int l_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);
            lfactor = l_ * l_ * (2.0 * l_ + 1.0) * (l_ + 1.0);

            std::visit(
               [&](auto&& p)
               {
                  spectrum = this->mValue.at(l_).transpose() *
                             p->profile(j, k).topRows(nN);
               },
               pInVarPolQ);

            this->mCmbSpectrum(l_) +=
               factor * lfactor * (spectrum.array().abs2()).real()(0);

            if (l_ == 1 and m_ == 0)
            {
               this->mAxialDipole = (l_ + 1.) * l_ * spectrum.real()(0);
            }
            if (l_ == 1 and m_ == 1)
            {
               this->mNonAxialDipole = (l_ + 1.) * l_ * spectrum(0);
            }
         }
      }
   }
   else
   {
      // Loop over harmonic degree l
      for (int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
      {
         int l_ = tRes.idx<Dimensions::Data::DAT3D>(k);
         lfactor = l_ * l_ * (2.0 * l_ + 1.0) * (l_ + 1.0);
         for (int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
         {
            int m_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);
            // m = 0, no factor of two
            if (m_ == 0)
            {
               factor = 1.0;
            }
            else
            {
               factor = 2.0;
            }

            std::visit(
               [&](auto&& p)
               {
                  spectrum = this->mValue.at(l_).transpose() *
                             p->profile(j, k).topRows(nN);
               },
               pInVarPolQ);

            this->mCmbSpectrum(l_) +=
               factor * lfactor * (spectrum.array().abs2()).real()(0);

            if (l_ == 1 and m_ == 0)
            {
               this->mAxialDipole = (l_ + 1.) * l_ * spectrum.real()(0);
            }
            if (l_ == 1 and m_ == 1)
            {
               this->mNonAxialDipole = (l_ + 1.) * l_ * spectrum(0);
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

void ShellDipolarityWriter::writeContent()
{
   // Create file
   this->preWrite();

// Get the "global" value from MPI code
#ifdef QUICC_MPI
   MPI_Allreduce(MPI_IN_PLACE, this->mCmbSpectrum.data(),
      this->mCmbSpectrum.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &this->mAxialDipole, 1, MPI_DOUBLE, MPI_SUM,
      MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &this->mNonAxialDipole, 1, MPI_DOUBLE_COMPLEX,
      MPI_SUM, MPI_COMM_WORLD);
#endif // QUICC_MPI

   using Tools::Formatter::ioFW;
   int ioPrec = 14;

   MHDFloat ro =
      this->mPhysical.find(NonDimensional::Upper1d::id())->second->value();
   this->mCmbSpectrum = this->mCmbSpectrum / ro;
   this->mAxialDipole = this->mAxialDipole / ro;
   this->mNonAxialDipole = this->mNonAxialDipole / ro;

   // Compute dipolarity
   this->mDipolarity =
      std::sqrt(this->mCmbSpectrum(1) / this->mCmbSpectrum.topRows(13).sum());

   // Check if the workflow allows IO to be performed
   if (QuICCEnv().allowsIO())
   {
      this->mFile << std::scientific;
      this->mFile << std::setprecision(ioPrec) << ioFW(ioPrec) << this->mTime
                  << "\t" << ioFW(ioPrec) << this->mDipolarity << "\t"
                  << this->mAxialDipole << "\t" << this->mNonAxialDipole.real()
                  << "\t" << this->mNonAxialDipole.imag();

      for (int count = 0; count < std::min(this->mCmbNl, static_cast<int>(this->mCmbSpectrum.size())); count++)
      {
         this->mFile << "\t" << this->mCmbSpectrum(count);
      }
      this->mFile << std::endl;
   }

   // Close file
   this->postWrite();

   // Abort if is NaN
   if (std::isnan(this->mCmbSpectrum.sum()))
   {
      QuICCEnv().abort("Shell dipolarity is NaN!");
   }
}

void ShellDipolarityWriter::setCmbTruncation(const int nL)
{
   if(nL < 0)
   {
      throw std::logic_error("Cannot set a negative truncation");
   }

   this->mCmbNl = nL;
}

void ShellDipolarityWriter::resetEnergy()
{
   this->mCmbSpectrum.setZero();
   this->mAxialDipole = 0.0;
   this->mNonAxialDipole = 0.0;
}

} // namespace Variable
} // namespace Io
} // namespace QuICC
