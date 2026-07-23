/**
 * @file SphereLuminosityWriter.cpp
 * @brief Source of the implementation of the ASCII Luminosity in a sphere
 */

 // Format of the  .dat file produced here:
// time,    total luminosity at ro,     nusselt at ro,
// where:
// total luminosity  = (backgraound luminosity) + (convective luminosity) (definition of Jones et al., 2011)
// nusselt           = (total luminosity) / (background luminosity)

// TODO: will probably need to come up with another definition


// System includes
//
#include <iomanip>
#include <stdexcept>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Io/Variable/SphereLuminosityWriter.hpp"
#include "QuICC/Io/Variable/Tags/Luminosity.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/dWnl.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

SphereLuminosityWriter::SphereLuminosityWriter(const std::string& prefix,
   const std::string& type, std::vector<std::shared_ptr<QuICC::DenseSM::Worland::RadialTorPolFunction>> pF) :
    IVariableAsciiWriter(prefix + Tags::Luminosity::BASENAME,
       Tags::Luminosity::EXTENSION, prefix + Tags::Luminosity::HEADER, type,
       Tags::Luminosity::VERSION, Dimensions::Space::SPECTRAL, EXTEND),
    mHasMOrdering(false),
    mLuminosity(std::numeric_limits<MHDFloat>::quiet_NaN()),
    mNusselt(std::numeric_limits<MHDFloat>::quiet_NaN()),
    mSb(std::numeric_limits<MHDFloat>::quiet_NaN()),
    mBoundary(0, 0)
{
   mpRhoTempKappa = pF[0];
   mpD1Sc = pF[1];
}

void SphereLuminosityWriter::init()
{
   this->mHasMOrdering =
      this->res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123);

   int m0 = -1;
   int l0 = -1;

   if (this->res()
          .cpu()
          ->dim(Dimensions::Transform::SPECTRAL)
          ->dim<Dimensions::Data::DAT3D>() > 0)
   {
      if (this->mHasMOrdering)
      {
         m0 = this->res()
                 .cpu()
                 ->dim(Dimensions::Transform::SPECTRAL)
                 ->idx<Dimensions::Data::DAT3D>(0);
         l0 = this->res()
                 .cpu()
                 ->dim(Dimensions::Transform::SPECTRAL)
                 ->idx<Dimensions::Data::DAT2D>(0, 0);
      }
      else
      {
         l0 = this->res()
                 .cpu()
                 ->dim(Dimensions::Transform::SPECTRAL)
                 ->idx<Dimensions::Data::DAT3D>(0);
         m0 = this->res()
                 .cpu()
                 ->dim(Dimensions::Transform::SPECTRAL)
                 ->idx<Dimensions::Data::DAT2D>(0, 0);
      }
   }

   // Background state
   Array rbArr(1);
   rbArr(0) = 1.0;
   MHDFloat rtk = (this->mpRhoTempKappa->evaluateLP(rbArr,0,0).array())(0);
   MHDFloat bg = -rtk*(this->mpD1Sc->evaluateLP(rbArr,0,0).array())(0);
   this->mSb = (4.0*Math::PI) * bg;

   // Look for l = 0, m = 0 mode
   if (m0 == 0 && l0 == 0)
   {
      Internal::Array grid = Internal::Array::Ones(1);
      int nN = this->res().sim().dim(Dimensions::Simulation::SIM1D,
         Dimensions::Space::SPECTRAL);
      Matrix poly(grid.size(), nN);
      Internal::Matrix ipoly(grid.size(), nN);
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::dWnl dwnl;
      dwnl.compute<MHDFloat>(poly, nN, 0, grid, Internal::Array(), ev::Set());
      this->mBoundary = poly.transpose();
      this->mBoundary *= -(4.0*Math::PI)*rtk / std::sqrt(4.0 * Math::PI);
   }
   else
   {
      this->mBoundary.resize(0, 0);
   }

   IVariableAsciiWriter::init();
}

void SphereLuminosityWriter::writeContent()
{
   scalar_iterator_range sRange = this->scalarRange();
   assert(std::distance(sRange.first, sRange.second) == 1);

   if (this->mBoundary.size() > 0)
   {
      this->mLuminosity = std::visit(
         [&](auto&& p)
         {
            return ((this->mSb + (this->mBoundary.transpose() *
                                                p->dom(0).total().profile(0, 0))
                                                .array()))
               .abs()(0, 0);
         },
         sRange.first->second);

         this->mNusselt = std::visit(
         [&](auto&& p)
         {
            return ((this->mSb + (this->mBoundary.transpose() *
                                                p->dom(0).total().profile(0, 0))
                                                .array())/this->mSb)
               .abs()(0, 0);
         },
         sRange.first->second);
   }
   else
   {
      this->mLuminosity = 0.0;
      this->mNusselt = 0.0;
   }

   // Create file
   this->preWrite();

// Get the "global" Kinetic energy from MPI code
#ifdef QUICC_MPI
   MPI_Allreduce(MPI_IN_PLACE, &this->mLuminosity, 1, MPI_DOUBLE, MPI_SUM,
      MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &this->mNusselt, 1, MPI_DOUBLE, MPI_SUM,
      MPI_COMM_WORLD);
#endif // QUICC_MPI

   using Tools::Formatter::ioFW;
   int ioPrec = 14;

   // Check if the workflow allows IO to be performed
   if (QuICCEnv().allowsIO())
   {
      this->mFile << std::scientific;
      this->mFile << std::setprecision(ioPrec) << ioFW(ioPrec) << this->mTime
                  << "\t" << ioFW(ioPrec) << this->mLuminosity << "\t" << ioFW(ioPrec) << this->mNusselt  << std::endl;
   }

   // Close file
   this->postWrite();

   // Abort if kinetic energy is NaN
   if (std::isnan(this->mLuminosity))
   {
      QuICCEnv().abort("Sphere Nusselt is NaN!");
   }
}

} // namespace Variable
} // namespace Io
} // namespace QuICC
