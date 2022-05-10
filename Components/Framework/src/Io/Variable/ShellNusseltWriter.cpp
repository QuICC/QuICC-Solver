/**
 * @file ShellNusseltWriter.cpp
 * @brief Source of the implementation of the ASCII Nusselt number in a spherical shell
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
#include "QuICC/Io/Variable/ShellNusseltWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/NonDimensional/Upper1D.hpp"
#include "QuICC/NonDimensional/Lower1D.hpp"
#include "QuICC/NonDimensional/Heating.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Io/Variable/Tags/Nusselt.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ShellNusseltWriter::ShellNusseltWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + Tags::Nusselt::BASENAME, Tags::Nusselt::EXTENSION, prefix + Tags::Nusselt::HEADER, type, Tags::Nusselt::VERSION, Dimensions::Space::SPECTRAL, EXTEND), mHasMOrdering(false), mNusselt(2), mBackground(2), mBoundary(0,0)
   {
   }

   ShellNusseltWriter::~ShellNusseltWriter()
   {
   }

   void ShellNusseltWriter::init()
   {
      this->mHasMOrdering = this->res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123);

      int m0, l0;
      if(this->mHasMOrdering)
      {
         m0 = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0);
         l0 = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,0);
      } else
      {
         l0 = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0);
         m0 = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,0);
      }

      // Look for l = 0, m = 0 mode
      if(m0 == 0 && l0 == 0)
      {
         auto ro = this->mPhysical.find(NonDimensional::Upper1D::id())->second->value();
         auto ri = this->mPhysical.find(NonDimensional::Lower1D::id())->second->value();
         auto a = (ro - ri)/2.0;

         this->mBackground.resize(2);
         int flag = this->mPhysical.find(NonDimensional::Heating::id())->second->value();
         // Internal heating
         if(flag == 0)
         {
            this->mBackground(0) = -ro;
            this->mBackground(1) = -ri;
         }
         else if(flag == 1)
         {
            auto c = ri*ro;
            this->mBackground(0) = -c*std::pow(ro,-2);
            this->mBackground(1) = -c*std::pow(ri,-2);
         }
         else
         {
            throw std::logic_error("Unknown background profile for spherical shell Nusselt writer");
         }

         int nN = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
         this->mBoundary.resize(nN, 2);
         for(int i = 0; i < this->mBoundary.rows(); i++)
         {
            this->mBoundary(i,0) = (2.0/a)*i*i/std::sqrt(4.0*Math::PI);
            this->mBoundary(i,1) = std::pow(-1,i+1)*this->mBoundary(i,0);
         }
      }
      else
      {
         this->mBoundary.resize(0,0);
         this->mBackground.resize(0);
      }

      IVariableAsciiWriter::init();
   }

   void ShellNusseltWriter::writeContent()
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);

      if(this->mBackground.size() > 0)
      {
         this->mNusselt = std::visit([&](auto&& p)->Array{return (this->mBackground + this->mBoundary.transpose()*p->dom(0).total().profile(0,0).real()).array()/this->mBackground.array();}, sRange.first->second);
      } else
      {
         this->mNusselt.setZero();
      }

      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mNusselt.data(), this->mNusselt.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      using Tools::Formatter::ioFW;
      int ioPrec = 14;

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::scientific;
         this->mFile << std::setprecision(ioPrec) << ioFW(ioPrec) << this->mTime << "\t" << ioFW(ioPrec) << this->mNusselt(0) << "\t" << ioFW(ioPrec) << this->mNusselt(1) << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mNusselt.sum()))
      {
         QuICCEnv().abort(99);

         throw std::logic_error("Spherical shell Nusselt is NaN!");
      }
   }

}
}
}
