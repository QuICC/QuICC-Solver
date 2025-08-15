/**
 * @file ShellLuminosityWriter.cpp
 * @brief Source of the implementation of the ASCII Luminosity in a spherical shell
 */

// Format of the  .dat file produced here:
// time,    total luminosity at ri,     nusselt at ri,      total luminosity at ro,     nusselt at ro,
// where:
// total luminosity  = (backgraound luminosity) + (convective luminosity) (definition of Jones et al., 2011)
// nusselt           = (total luminosity) / (background luminosity)

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
#include "QuICC/Io/Variable/ShellLuminosityWriter.hpp"

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "Types/Math.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Heating.hpp"
#include "QuICC/NonDimensional/Beta.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Io/Variable/Tags/Luminosity.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ShellLuminosityWriter::ShellLuminosityWriter(const std::string& prefix, const std::string& type, std::vector<std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction>> pF)
      : IVariableAsciiWriter(prefix + Tags::Luminosity::BASENAME, Tags::Luminosity::EXTENSION, prefix + Tags::Luminosity::HEADER, type, Tags::Luminosity::VERSION, Dimensions::Space::SPECTRAL, EXTEND), mHasMOrdering(false), mNusselt(2), mLuminosity(2), mBackground(2), mBoundary(0,0)
   {
      assert(pF.size() >= 2 && "ShellLuminosityWriter requires at least 2 elements in pF vector");
      mpRhoTempKappa = pF[0];
      mpD1Sc = pF[1];
   }

   ShellLuminosityWriter::~ShellLuminosityWriter()
   {
   }

   void ShellLuminosityWriter::init()
   {
      this->mHasMOrdering = this->res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123);
      const auto& tRes = *this->res().cpu()->dim(Dimensions::Transform::SPECTRAL);

      int m0, l0;
      if(this->mHasMOrdering)
      {
         m0 = tRes.idx<Dimensions::Data::DAT3D>(0);
         l0 = tRes.idx<Dimensions::Data::DAT2D>(0,0);
      } else
      {
         l0 = tRes.idx<Dimensions::Data::DAT3D>(0);
         m0 = tRes.idx<Dimensions::Data::DAT2D>(0,0);
      }

      // Look for l = 0, m = 0 mode
      if(m0 == 0 && l0 == 0)
      {
         auto ro = this->mPhysical.find(NonDimensional::Upper1d::id())->second->value();
         auto ri = this->mPhysical.find(NonDimensional::Lower1d::id())->second->value();
         auto a = (ro - ri)/2.0;

         Internal::Array rbArr(2);
         rbArr(0) = ro;
         rbArr(1) = ri;


         this->mBackground.resize(2);
         int flag = this->mPhysical.find(NonDimensional::Heating::id())->second->value();
         // Internal heating
         if(flag == 0)
         {
            Internal::Array bgArray = -(this->mpRhoTempKappa->evaluate(rbArr,0,0).array()) * (this->mpD1Sc->evaluate(rbArr,0,0).array());
            //this->mBackground(0) = -ro;
            //this->mBackground(1) = -ri;
            this->mBackground(0) = bgArray(0) * (4.0*Math::PI)*ro*ro;
            this->mBackground(1) = bgArray(1) * (4.0*Math::PI)*ri*ri;    
         }
         else if(flag == 1)
         {
            throw std::logic_error("Unknown background profile for spherical shell Luminosity writer. Potentially outdated option");
         }
         else if(flag == 2 || flag ==3)
         {            
            throw std::logic_error("Unknown background profile for spherical shell Luminosity writer. Potentially outdated option");
         }
         else
         {
            throw std::logic_error("Unknown background profile for spherical shell Luminosity writer");
         }

         int nN = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
         this->mBoundary.resize(nN, 2);
         for(int i = 0; i < this->mBoundary.rows(); i++)
         {
            //this->mBoundary(i,0) = (2.0/a)*i*i/std::sqrt(4.0*Math::PI);
            //this->mBoundary(i,1) = std::pow(-1,i+1)*this->mBoundary(i,0);
            this->mBoundary(i,0) = (2.0/a)*i*i/std::sqrt(4.0*Math::PI)   * (4.0*Math::PI)*ro*ro*(-this->mpRhoTempKappa->evaluate(rbArr,0,0).array()(0));
            this->mBoundary(i,1) = std::pow(-1,i+1)*this->mBoundary(i,0) * (4.0*Math::PI)*ri*ri*(-this->mpRhoTempKappa->evaluate(rbArr,0,0).array()(1)); 
         }
      }
      else
      {
         this->mBoundary.resize(0,0);
         this->mBackground.resize(0);
      }

      IVariableAsciiWriter::init();
   }

   void ShellLuminosityWriter::writeContent()
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);

      if(this->mBackground.size() > 0)
      {
         this->mNusselt = std::visit([&](auto&& p)->Array{return (this->mBackground + this->mBoundary.transpose()*p->dom(0).total().profile(0,0).real()).array()/this->mBackground.array();}, sRange.first->second);
         this->mLuminosity = std::visit([&](auto&& p)->Array{return (this->mBackground + this->mBoundary.transpose()*p->dom(0).total().profile(0,0).real()).array();}, sRange.first->second);
      } else
      {
         this->mNusselt.setZero();
         this->mLuminosity.setZero();
      }

      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mNusselt.data(), this->mNusselt.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, this->mLuminosity.data(), this->mLuminosity.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      using Tools::Formatter::ioFW;
      int ioPrec = 14;

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::scientific;
         //this->mFile << std::setprecision(ioPrec) << ioFW(ioPrec) << this->mTime << "\t" << ioFW(ioPrec) << this->mNusselt(0) << "\t" << ioFW(ioPrec) << this->mNusselt(1) << std::endl;
         this->mFile << std::setprecision(ioPrec) << ioFW(ioPrec) << this->mTime << "\t" << ioFW(ioPrec) << this->mLuminosity(0) << "\t" << ioFW(ioPrec) << this->mNusselt(0) << "\t" << ioFW(ioPrec) << this->mLuminosity(1) << "\t" << ioFW(ioPrec) << this->mNusselt(1) << std::endl;

      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mNusselt.sum()) || std::isnan(this->mLuminosity.sum()))
      {
         QuICCEnv().abort("Spherical shell Luminosity or Nusselt is NaN!");
      }
   }

}
}
}
