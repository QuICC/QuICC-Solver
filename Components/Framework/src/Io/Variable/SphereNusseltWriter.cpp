/**
 * @file SphereNusseltWriter.cpp
 * @brief Source of the implementation of the ASCII Nusselt number in a sphere
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
#include "QuICC/Io/Variable/SphereNusseltWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Io/Variable/Tags/Nusselt.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   SphereNusseltWriter::SphereNusseltWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + Tags::Nusselt::BASENAME, Tags::Nusselt::EXTENSION, prefix + Tags::Nusselt::HEADER, type, Tags::Nusselt::VERSION, Dimensions::Space::SPECTRAL, EXTEND), mHasMOrdering(false), mNusselt(std::numeric_limits<MHDFloat>::quiet_NaN()), mTb(std::numeric_limits<MHDFloat>::quiet_NaN()), mOrigin(0,0)
   {
   }

   SphereNusseltWriter::~SphereNusseltWriter()
   {
   }

   void SphereNusseltWriter::init()
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

      // Background state
      this->mTb = 0.5;

      // Look for l = 0, m = 0 mode
      if(m0 == 0 && l0 == 0)
      {
         internal::Array grid = internal::Array::Zero(1);
         int nN = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
         Matrix poly(grid.size(), nN);
         internal::Matrix  ipoly(grid.size(), nN);
         namespace ev = Polynomial::Worland::Evaluator;
         Polynomial::Worland::Wnl wnl;
         wnl.compute<MHDFloat>(poly, nN, 0, grid, internal::Array(), ev::Set());
         this->mOrigin = poly.transpose();
         this->mOrigin /= std::sqrt(4.0*Math::PI);
      } else
      {
         this->mOrigin.resize(0,0);
      }

      IVariableAsciiWriter::init();
   }

   void SphereNusseltWriter::writeContent()
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);

      if(this->mOrigin.size() > 0)
      {
         this->mNusselt = std::visit([&](auto&& p){return (this->mTb/(this->mTb + (this->mOrigin.transpose()*p->dom(0).total().profile(0,0)).array())).abs()(0,0);}, sRange.first->second);
      } else
      {
         this->mNusselt = 0.0;
      }

      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &this->mNusselt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      using Tools::Formatter::ioFW;
      int ioPrec = 14;

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::scientific;
         this->mFile << std::setprecision(ioPrec) << ioFW(ioPrec) << this->mTime << "\t" << ioFW(ioPrec) << this->mNusselt << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mNusselt))
      {
         QuICCEnv().abort("Sphere Nusselt is NaN!");
      }
   }

}
}
}
