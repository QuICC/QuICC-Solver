/**
 * @file CartesianCflWrapper.cpp
 * @brief Source of the CFL constraint wrapper in a Cartesian geometry
 */

// Debug includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Diagnostics/CartesianCflWrapper.hpp"

// Project includes
//

namespace QuICC {

namespace Diagnostics {

   CartesianCflWrapper::CartesianCflWrapper(const SharedIVectorWrapper spVelocity)
      : ICflWrapper(spVelocity), mcCourant(0.65)
   {
   }

   CartesianCflWrapper::~CartesianCflWrapper()
   {
   }

   void CartesianCflWrapper::init(const std::vector<Array>& mesh)
   {
      // Initialize the mesh
      this->initMesh(mesh);
   }

   void CartesianCflWrapper::initMesh(const std::vector<Array>& mesh)
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

   Matrix CartesianCflWrapper::initialCfl() const
   {
      Matrix cfl = this->cfl();

      const Array& dx1 = this->mMeshSpacings.at(0);
      const Array& dx2 = this->mMeshSpacings.at(1);
      const Array& dx3 = this->mMeshSpacings.at(2);

      // Assume a velocity of 100 to avoid problems with "zero" starting values
      MHDFloat newCfl;
      int idx;
      newCfl = this->mcCourant*dx1.minCoeff(&idx)/100.;
      if(newCfl < cfl(0,1))
      {
         cfl(0,1) = newCfl;
         cfl(1,1) = dx1(idx);
      }
      newCfl = this->mcCourant*dx2.minCoeff(&idx)/100.;
      if(newCfl < cfl(0,2))
      {
         cfl(0,2) = newCfl;
         cfl(1,2) = dx2(idx);
      }
      newCfl = this->mcCourant*dx3.minCoeff(&idx)/100.;
      if(newCfl < cfl(0,3))
      {
         cfl(0,3) = newCfl;
         cfl(1,3) = dx3(idx);
      }

      this->updateCflMatrix(cfl);

      return cfl;
   }

   Matrix CartesianCflWrapper::cfl() const
   {
      // Compute most stringent CFL condition
      MHDFloat newCfl;
      Matrix cfl = Matrix::Constant(2,4, std::numeric_limits<MHDFloat>::max());

      const Array& dx1 = this->mMeshSpacings.at(0);
      const Array& dx2 = this->mMeshSpacings.at(1);
      const Array& dx3 = this->mMeshSpacings.at(2);

      int nK = this->mspVelocity->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();

      // CFL from first component
      for(int k = 0; k < nK; ++k)
      {
         int k_ = this->mspVelocity->res().cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(k);
         newCfl = dx1(k_)/this->mspVelocity->one().slice(k).array().abs().maxCoeff();
         if(newCfl < cfl(0,1))
         {
            cfl(0,1) = newCfl;
            cfl(1,1) = 0;
         }
      }

      // CFL from second component
      for(int k = 0; k < nK; ++k)
      {
         int nJ = this->mspVelocity->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(k);
         for(int j = 0; j < nJ; ++j)
         {
            int j_ = this->mspVelocity->res().cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(j,k);
            newCfl = dx2(j_)/this->mspVelocity->two().profile(j,k).array().abs().maxCoeff();
            if(newCfl < cfl(0,2))
            {
               cfl(0,2) = newCfl;
               cfl(1,2) = 0;
            }
         }
      }

      // CFL from second component
      for(int k = 0; k < nK; ++k)
      {
         int nJ = this->mspVelocity->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(k);
         for(int j = 0; j < nJ; ++j)
         {
            int nI = this->mspVelocity->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DATF1D>(j,k);
            for(int i = 0; i < nI; ++i)
            {
               newCfl = dx3(i)/std::abs(this->mspVelocity->three().point(i,j,k));
               if(newCfl < cfl(0,3))
               {
                  cfl(0,3) = newCfl;
                  cfl(1,3) = 0;
               }
            }
         }
      }

      cfl.row(0).tail(cfl.cols()-1).array() *= this->mcCourant;
      this->updateCflMatrix(cfl);

      return cfl;
   }

}
}
