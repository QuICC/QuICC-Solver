/**
 * @file CartesianCfl.cpp
 * @brief Source of the CFL constraint wrapper in a Cartesian geometry
 */

// System includes
//

// Project includes
//
#include "QuICC/Diagnostics/CartesianCfl.hpp"

namespace QuICC {

namespace Diagnostics {

   CartesianCfl::CartesianCfl(const std::map<std::size_t,NonDimensional::SharedINumber>& params, const MHDFloat courant)
      : ICflWrapper(courant), mVelId(0)
   {
      this->mIsActive = true;
   }

   void CartesianCfl::defineVelocity(const std::size_t velId)
   {
      this->mVelId = velId;
      this->mFields.try_emplace(this->mVelId, nullptr);
   }

   void CartesianCfl::init(const std::vector<Array>& mesh)
   {
      // Initialize the mesh
      this->initMesh(mesh);
   }

   void CartesianCfl::initMesh(const std::vector<Array>& mesh)
   {
      // Compute the mesh spacings
      this->mMeshSpacings.reserve(mesh.size());
      // Loop over all dimensions
      for(size_t i = 0; i < mesh.size(); i++)
      {
         // Create storage
         this->mMeshSpacings.push_back(Array(mesh.at(i).size()));

         // Extract minimal spacing for each grid in current direction
         for(int j = 0; j < mesh.at(i).size(); ++j)
         {
            // Get internal points grid spacing
            if(j > 0 && j < mesh.at(i).size() - 1)
            {
               this->mMeshSpacings.back()(j) = std::min(std::abs(mesh.at(i)(j) - mesh.at(i)(j-1)), std::abs(mesh.at(i)(j) - mesh.at(i)(j+1)));

            // Get left endpoint grid spacing
            } else if(j > 0)
            {
               this->mMeshSpacings.back()(j) = std::abs(mesh.at(i)(j) - mesh.at(i)(j-1));

            // Get right endpoint grid spacing
            } else
            {
               this->mMeshSpacings.back()(j) = std::abs(mesh.at(i)(j) - mesh.at(i)(j+1));
            }
         }
      }
   }

   Matrix CartesianCfl::initialCfl() const
   {
      Matrix cfl = this->cfl();

      const Array& dx1 = this->mMeshSpacings.at(0);
      const Array& dx2 = this->mMeshSpacings.at(1);
      const Array& dx3 = this->mMeshSpacings.at(2);

      // Assume a velocity of 100 to avoid problems with "zero" starting values
      MHDFloat newCfl;
      int idx;
      newCfl = this->mcCourant*dx1.minCoeff(&idx)/100.;
      if(newCfl < cfl(0,0))
      {
         cfl(0,0) = newCfl;
         cfl(1,0) = dx1(idx);
      }
      newCfl = this->mcCourant*dx2.minCoeff(&idx)/100.;
      if(newCfl < cfl(0,1))
      {
         cfl(0,1) = newCfl;
         cfl(1,1) = dx2(idx);
      }
      newCfl = this->mcCourant*dx3.minCoeff(&idx)/100.;
      if(newCfl < cfl(0,2))
      {
         cfl(0,2) = newCfl;
         cfl(1,2) = dx3(idx);
      }

      return cfl;
   }

   Matrix CartesianCfl::cfl() const
   {
      // Compute most stringent CFL condition
      MHDFloat newCfl;
      Matrix cfl = Matrix::Constant(2,3, std::numeric_limits<MHDFloat>::max());

      const Array& dx1 = this->mMeshSpacings.at(0);
      const Array& dx2 = this->mMeshSpacings.at(1);
      const Array& dx3 = this->mMeshSpacings.at(2);

      const auto& vel = this->mFields.at(this->mVelId);

      int nK = vel->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();

      // CFL from first component
      for(int k = 0; k < nK; ++k)
      {
         int k_ = vel->res().cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(k);
         newCfl = dx1(k_)/vel->one().slice(k).array().abs().maxCoeff();
         if(newCfl < cfl(0,0))
         {
            cfl(0,0) = newCfl;
            cfl(1,0) = 0;
         }
      }

      // CFL from second component
      for(int k = 0; k < nK; ++k)
      {
         int nJ = vel->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(k);
         for(int j = 0; j < nJ; ++j)
         {
            int j_ = vel->res().cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(j,k);
            newCfl = dx2(j_)/vel->two().profile(j,k).array().abs().maxCoeff();
            if(newCfl < cfl(0,1))
            {
               cfl(0,1) = newCfl;
               cfl(1,1) = 0;
            }
         }
      }

      // CFL from second component
      for(int k = 0; k < nK; ++k)
      {
         int nJ = vel->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(k);
         for(int j = 0; j < nJ; ++j)
         {
            int nI = vel->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DATF1D>(j,k);
            for(int i = 0; i < nI; ++i)
            {
               newCfl = dx3(i)/std::abs(vel->three().point(i,j,k));
               if(newCfl < cfl(0,2))
               {
                  cfl(0,2) = newCfl;
                  cfl(1,2) = 0;
               }
            }
         }
      }

      cfl.row(0).array() *= this->mcCourant;

      return cfl;
   }

}
}
