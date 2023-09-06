/** 
 * @file IRegular2DBuilder.cpp
 * @brief Source of the generic regulard 2D scheme implementation
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SpatialScheme/2D/IRegular2DBuilder.hpp"
#include "QuICC/SpatialScheme/Tools/Uniform3D1D.hpp"
#include "QuICC/SpatialScheme/Tools/Uniform3D2D.hpp"
#include "QuICC/Transform/Setup/Uniform.hpp"

namespace QuICC {

namespace SpatialScheme {

   IRegular2DBuilder::IRegular2DBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options)
      : IBuilder(dim.size(), purpose, options), mI(dim(0)), mJ(dim(1))
   {
      assert(dim.size() == 3);
   }

   ArrayI IRegular2DBuilder::resolution() const
   {
      ArrayI space(this->dims());
      space << this->mI, this->mJ;

      return space;
   }

   int IRegular2DBuilder::fillIndexes(const Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const std::vector<int>& id, const std::vector<int>& bins)
   {
      throw std::logic_error("2D scheme splitting not fully implemented");

      // Safety assertions
      assert( id.size() > 0 );
      assert( bins.size() > 0 );
      assert( id.size() == bins.size() );

      // Make sure we start with empty indexes
      fwd1D.clear();
      bwd1D.clear();
      idx2D.clear();

      // Multimap for the modes
      std::multimap<int,int> modes;

      this->split(modes, id, bins, transId);

      // Truncation tools
      auto spTools = this->truncationTools(transId);

      // Fill indexes for 2D and 3D
      spTools->fillIndexes2D3D(idx2D, idx3D, modes);

      // Fill indexes for 1D
      spTools->fillIndexes1D(fwd1D, bwd1D, idx3D, this->dim(transId, Dimensions::Data::DATF1D), this->dim(transId, Dimensions::Data::DATB1D));

      // Set status (0 for success, 1 for failure)
      int status = 0;
      if(modes.size() == 0)
      {
         status = 1;
      }

      return status;
   }

   int IRegular2DBuilder::splittableTotal(const  Dimensions::Transform::Id transId, Splitting::Locations::Id flag)
   {
      if(flag == Splitting::Locations::FIRST || flag == Splitting::Locations::COUPLED2D)
      {
         return this->dim(transId, Dimensions::Data::DAT2D);
      }
      
      // If condition has not been mached
      throw std::logic_error("Can't split in any other dimension than first dimension in 2D regular case");

      return -1;
   }

   std::shared_ptr<Tools::IBase> IRegular2DBuilder::truncationTools(const Dimensions::Transform::Id transId) const
   {
      // Setup truncation tools
      std::shared_ptr<Tools::IBase> spTools;

      if(transId == Dimensions::Transform::TRA1D || transId == Dimensions::Transform::SPECTRAL)
      {
         spTools = std::make_shared<Tools::Uniform3D1D>();
      }
      else if(transId == Dimensions::Transform::TRA2D)
      {
         spTools = std::make_shared<Tools::Uniform3D2D>();
      }
      else if(transId == Dimensions::Transform::TRA3D)
      {
         throw std::logic_error("Tried to work on third dimension in 2D case");
      }

      return spTools;
   }

   void IRegular2DBuilder::split(std::multimap<int,int>& modes, const std::vector<int>& id, const std::vector<int>& bins, const Dimensions::Transform::Id transId)
   {
      auto spTools = this->truncationTools(transId);
      std::vector<int> n0_, nN_;

      auto n2D = this->dim(transId, Dimensions::Data::DAT2D);
      auto n3D = this->dim(transId, Dimensions::Data::DAT3D);

      // Create list of modes
      spTools->buildBalancedMap(modes, n2D, n3D, id, bins, n0_, nN_);
   }

   void IRegular2DBuilder::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);

      // Set second transform cost
      this->setCost(1.0, Dimensions::Transform::TRA2D);
   }

   void IRegular2DBuilder::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);

      // Set second transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA2D);
   }

   void IRegular2DBuilder::setMemoryScore()
   {
      // Set first transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA1D);

      // Set second transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA2D);
   }

   void IRegular2DBuilder::setDimensions()
   {
      //
      // Set transform space sizes
      //
      ArrayI traSize(2);
      traSize(0) = this->mesher().nSpec1D();
      traSize(1) = this->mesher().nSpec2D();
      this->setTransformSpace(traSize);

      //
      // Initialise spectral space
      //

      // Initialise forward dimension of first transform
      this->setDimension(this->mesher().nSpec1D(), Dimensions::Transform::SPECTRAL, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(this->mesher().nSpec1D(), Dimensions::Transform::SPECTRAL, Dimensions::Data::DATB1D);

      // Initialise second dimension of first transform
      this->setDimension(this->mesher().nSpec2D(), Dimensions::Transform::SPECTRAL, Dimensions::Data::DAT2D);

      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->setDimension(this->mesher().nPhys1D(), Dimensions::Transform::TRA1D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(this->mesher().nDealias1D(), Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D);

      // Initialise second dimension of first transform
      this->setDimension(this->mesher().nSpec2D(), Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

      //
      // Initialise second transform
      //

      // Initialise forward dimension of second transform
      this->setDimension(this->mesher().nPhys2D(), Dimensions::Transform::TRA2D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of second transform
      this->setDimension(this->mesher().nDealias2D(), Dimensions::Transform::TRA2D, Dimensions::Data::DATB1D);

      // Initialise second dimension of second transform
      this->setDimension(this->mesher().nPhys1D(), Dimensions::Transform::TRA2D, Dimensions::Data::DAT2D);
   }

} // SpatialScheme
} // QuICC
