/**
 * @file SphereAngularMomentumWriter.cpp
 * @brief Source of the implementation of the ASCII angular momentum number in a sphere
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
#include "QuICC/Io/Variable/SphereAngularMomentumWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Io/Variable/Tags/AngularMomentum.hpp"
#include "QuICC/Polynomial/Worland/Operators.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   SphereAngularMomentumWriter::SphereAngularMomentumWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + Tags::AngularMomentum::BASENAME, Tags::AngularMomentum::EXTENSION, prefix + Tags::AngularMomentum::HEADER, type, Tags::AngularMomentum::VERSION, Dimensions::Space::SPECTRAL, EXTEND), mHasMOrdering(false), mHasM0(false), mHasM1(false), mM0j(-1), mM0k(-1), mM1j(-1), mM1k(-1), mMomentum(3)
   {
   }

   SphereAngularMomentumWriter::~SphereAngularMomentumWriter()
   {
   }

   void SphereAngularMomentumWriter::init()
   {
      this->mHasMOrdering = this->res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123);

      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int m_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

               if(l_ == 1 && m_ == 0)
               {
                  this->mM0j = j;
                  this->mM0k = k;
                  this->mHasM0 = true;
               } else if(l_ == 1 && m_ == 1)
               {
                  this->mM1j = j;
                  this->mM1k = k;
                  this->mHasM1 = true;
               }
            }
         }
      } else
      {
         // Loop over harmonic degree l
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int m_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
               if(l_ == 1 && m_ == 0)
               {
                  this->mM0j = j;
                  this->mM0k = k;
                  this->mHasM0 = true;
               } else if(l_ == 1 && m_ == 1)
               {
                  this->mM1j = j;
                  this->mM1k = k;
                  this->mHasM1 = true;
               }
            }
         }
      }

      // Compute operator if required
      if(this->mHasM0 || this->mHasM1)
      {
         int nN = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
         internal::Matrix iop;
         Polynomial::Worland::Operators::integrateRpWnl(iop, 1, 3, nN);
         this->mOp = iop.cast<MHDFloat>();
         assert(this->mOp.rows() == nN && this->mOp.cols() == 1);
      }

      IVariableAsciiWriter::init();
   }

   void SphereAngularMomentumWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      // get iterator to field
      vector_iterator vIt;
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);

      auto& field = vRange.first->second;
      assert(std::visit([&](auto&& p)->bool{return (p->dom(0).res().sim().ss().spectral().ONE() == FieldComponents::Spectral::TOR);}, field));
      assert(std::visit([&](auto&& p)->bool{return (p->dom(0).res().sim().ss().spectral().TWO() == FieldComponents::Spectral::POL);}, field));

      ArrayZ mom;
      this->mMomentum.setZero();
      if(this->mHasM1)
      {
         MHDFloat c = std::sqrt(32.0*Math::PI/3.0);
         std::visit([&](auto&& f){
            mom = (this->mOp.transpose()*f->dom(0).total().comp(FieldComponents::Spectral::TOR).profile(this->mM1j, this->mM1k));
         }, field);
         this->mMomentum(0) = -c*mom(0).real();
         this->mMomentum(1) = c*mom(0).imag();
      }

      if(this->mHasM0)
      {
         MHDFloat c = std::sqrt(16.0*Math::PI/3.0);
         std::visit([&](auto&& f){
            mom = (this->mOp.transpose()*f->dom(0).total().comp(FieldComponents::Spectral::TOR).profile(this->mM0j, this->mM0k));
         }, field);
         this->mMomentum(2) = c*mom(0).real();
      }
   }

   void SphereAngularMomentumWriter::writeContent()
   {
      // Create file
      this->preWrite();

      // Get the "global" value from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mMomentum.data(), this->mMomentum.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      using Tools::Formatter::ioFW;
      int ioPrec = 14;

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::scientific;
         this->mFile << std::setprecision(ioPrec) << ioFW(ioPrec) << this->mTime << "\t" << ioFW(ioPrec) << this->mMomentum.norm();
         for(int i = 0; i < this->mMomentum.size(); i++)
         {
            this->mFile << "\t" << ioFW(ioPrec) << this->mMomentum(i);
         }
         this->mFile << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if is NaN
      if(std::isnan(this->mMomentum.sum()))
      {
         QuICCEnv().abort("Sphere angular momentum is NaN!");
      }
   }

}
}
}
