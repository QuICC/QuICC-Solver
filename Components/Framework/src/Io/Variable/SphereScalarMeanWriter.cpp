/**
 * @file SphereScalarMeanWriter.cpp
 * @brief Source of the implementation of the scalar field mean value writer in a sphere
 */

// System includes
//
#include <iomanip>
#include <stdexcept>

// Project includes
//
#include "QuICC/Io/Variable/SphereScalarMeanWriter.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Io/Variable/Tags/Mean.hpp"
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   SphereScalarMeanWriter::SphereScalarMeanWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + Tags::Mean::BASENAME, Tags::Mean::EXTENSION, prefix + Tags::Mean::HEADER, type, Tags::Mean::VERSION, Dimensions::Space::SPECTRAL, EXTEND), mHasMOrdering(false), mMean(std::numeric_limits<MHDFloat>::quiet_NaN()), mVolume(std::numeric_limits<MHDFloat>::quiet_NaN()), mProjector(0,0)
   {
   }

   void SphereScalarMeanWriter::init()
   {
      this->mHasMOrdering = this->res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123);

      int m0 = -1;
      int l0 = -1;

      if(this->res().cpu()->dim(Dimensions::Transform::SPECTRAL)->dim<Dimensions::Data::DAT3D>() > 0)
      {
         if(this->mHasMOrdering)
         {
            m0 = this->res().cpu()->dim(Dimensions::Transform::SPECTRAL)->idx<Dimensions::Data::DAT3D>(0);
            l0 = this->res().cpu()->dim(Dimensions::Transform::SPECTRAL)->idx<Dimensions::Data::DAT2D>(0,0);
         } else
         {
            l0 = this->res().cpu()->dim(Dimensions::Transform::SPECTRAL)->idx<Dimensions::Data::DAT3D>(0);
            m0 = this->res().cpu()->dim(Dimensions::Transform::SPECTRAL)->idx<Dimensions::Data::DAT2D>(0,0);
         }
      }

      // Normalize by sphere volume: 4/3*pi*r_o^3
      this->mVolume = (4.0/3.0)*Math::PI;
      // Normalize by horizontal integral / spherical harmonic Y_0^0 normalization
      this->mVolume /= (4.0*Math::PI)/(2.0*std::sqrt(Math::PI));

      // Look for l = 0, m = 0 mode
      if(m0 == 0 && l0 == 0)
      {
         int nN = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
         int nR = 2*nN;
         internal::Array igrid,iweights;
         Polynomial::Quadrature::LegendreRule lrule;
         lrule.computeQuadrature(igrid, iweights, nR);
         igrid = (igrid.array() + MHD_MP(1.0))/MHD_MP(2.0);
         iweights *= MHD_MP(0.5);
         Matrix poly(igrid.size(), nN);
         internal::Matrix  ipoly(igrid.size(), nN);
         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;
         wnl.compute<internal::MHDFloat>(ipoly, nN, 0, igrid, (iweights.array()*igrid.array().pow(2)).matrix(), ev::Set());
         this->mProjector = ipoly.cast<MHDFloat>().transpose();
         this->mProjector /= this->mVolume;

         this->mBg = ArrayZ::Zero(nN);
         // Temperature background state T_b = 1/2(1-r^2)
         if(this->scalarRange().first->first == PhysicalNames::Temperature::id())
         {
            this->mBg(0) = Math::PI/std::sqrt(8.0);
            this->mBg(1) = -Math::PI/4.0;
         }
      }
      else
      {
         this->mProjector.resize(0,0);
      }

      IVariableAsciiWriter::init();
   }

   void SphereScalarMeanWriter::writeContent()
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);

      if(this->mProjector.size() > 0)
      {
         this->mMean = std::visit(
               [&](auto&& p)
               {
                  MHDComplex val = (this->mProjector.transpose()*(p->dom(0).total().profile(0,0)+this->mBg)).array().sum();
                  return val.real();
               }, sRange.first->second);
      }
      else
      {
         this->mMean = 0.0;
      }

      // Create file
      this->preWrite();

      // Get the "global" mean value from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &this->mMean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      using Tools::Formatter::ioFW;
      int ioPrec = 14;

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::scientific;
         this->mFile << std::setprecision(ioPrec) << ioFW(ioPrec) << this->mTime << "\t" << ioFW(ioPrec) << this->mMean << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if mean value is NaN
      if(std::isnan(this->mMean))
      {
         QuICCEnv().abort("Sphere scalar mean is NaN!");
      }
   }

}
}
}
