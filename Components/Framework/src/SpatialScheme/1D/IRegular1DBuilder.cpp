/** 
 * @file IRegular1DBuilder.cpp
 * @brief Source of the generic regular 1D scheme implementation
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SpatialScheme/1D/IRegular1DBuilder.hpp"
#include "QuICC/SpatialScheme/Tools/Uniform3D1D.hpp"
#include "QuICC/Transform/Setup/Uniform.hpp"

namespace QuICC {

namespace SpatialScheme {

   IRegular1DBuilder::IRegular1DBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options)
      : IBuilder(dim.size(), purpose, options), mI(dim(0))
   {
      assert(dim.size() == 3);
   }

   ArrayI IRegular1DBuilder::resolution() const
   {
      ArrayI space(this->dims());
      space << this->mI;

      return space;
   }

   int IRegular1DBuilder::fillIndexes(Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const std::vector<int>& id, const std::vector<int>& bins)
   {
      // Safety assertions
      assert( id.size() > 0 );
      assert( bins.size() > 0 );
      assert( id.size() == bins.size() );

      // Assert for right transform (1D case)
      assert(transId == Dimensions::Transform::TRA1D);

      // Set unused third dimension
      idx3D.resize(1);
      idx3D.setZero();

      // Clear second dimension
      idx2D.clear();
      idx2D.push_back(idx3D);

      // Make sure we start with empty indexes
      fwd1D.clear();
      bwd1D.clear();

      // Create single forward storage for indexes
      fwd1D.push_back(ArrayI(this->dim(transId, Dimensions::Data::DATF1D)));

      // Fill array with indexes
      for(int i = 0; i < fwd1D.at(0).size(); i++)
      {
         fwd1D.at(0)(i) = i;
      }

      // Create single backward storage for indexes
      bwd1D.push_back(ArrayI(this->dim(transId, Dimensions::Data::DATB1D)));

      // Fill array with indexes
      for(int i = 0; i < bwd1D.at(0).size(); i++)
      {
         bwd1D.at(0)(i) = i;
      }

      return 0;
   }

   int IRegular1DBuilder::splittableTotal(Dimensions::Transform::Id transId, Splitting::Locations::Id flag)
   {
      throw std::logic_error("There is no splitting algorithm for 1D problems!");
   }

   std::shared_ptr<Tools::IBase> IRegular1DBuilder::truncationTools(const Dimensions::Transform::Id transId) const
   {
      // Setup truncation tools
      std::shared_ptr<Tools::IBase> spTools;

      if(transId == Dimensions::Transform::TRA1D || transId == Dimensions::Transform::SPECTRAL)
      {
         spTools = std::make_shared<Tools::Uniform3D1D>();
      }
      else if(transId == Dimensions::Transform::TRA2D)
      {
         throw std::logic_error("Tried to work on second dimension in 1D case");
      }
      else if(transId == Dimensions::Transform::TRA3D)
      {
         throw std::logic_error("Tried to work on third dimension in 1D case");
      }

      return spTools;
   }

   void IRegular1DBuilder::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);
   }

   void IRegular1DBuilder::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);
   }

   void IRegular1DBuilder::setMemoryScore()
   {
      // Set first transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA1D);
   }

   void IRegular1DBuilder::setDimensions()
   {
      //
      // Set transform space sizes
      //
      ArrayI traSize(1);
      traSize(0) = this->mesher().nSpec1D();
      this->setTransformSpace(traSize);

      //
      // Initialise spectral space
      //

      // Initialise forward dimension of first transform
      this->setDimension(this->mesher().nSpec1D(), Dimensions::Transform::SPECTRAL, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(this->mesher().nSpec1D(), Dimensions::Transform::SPECTRAL, Dimensions::Data::DATB1D);

      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->setDimension(this->mesher().nPhys1D(), Dimensions::Transform::TRA1D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(this->mesher().nDealias1D(), Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D);
   }

} // SpatialScheme
} // QuICC
