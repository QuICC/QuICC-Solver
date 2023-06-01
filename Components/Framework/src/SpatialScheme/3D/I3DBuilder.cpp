/** 
 * @file I3DBuilder.cpp
 * @brief Source of the regular 3D scheme builder
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SpatialScheme/3D/I3DBuilder.hpp"
#include "QuICC/SpatialScheme/Tools/IBase.hpp"
#include "QuICC/SpatialScheme/Tools/Uniform3D1D.hpp"
#include "QuICC/SpatialScheme/Tools/Uniform3D2D.hpp"
#include "QuICC/SpatialScheme/Tools/Uniform3D3D.hpp"
#include "QuICC/Transform/Setup/Uniform.hpp"

namespace QuICC {

namespace SpatialScheme {

   I3DBuilder::I3DBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options)
      : IBuilder(dim.size(), purpose, options)
   {
      assert(dim.size() == 3);
   }

   int I3DBuilder::checkModes(const std::multimap<int,int>& modes, const Dimensions::Transform::Id transId) const
   {
      int status = 0;
      if(modes.size() == 0)
      {
         status = 1;
      }

      return status;
   }

   int I3DBuilder::fillIndexes(const Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const std::vector<int>& id, const std::vector<int>& bins)
   {
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
      int status = checkModes(modes, transId);
      return status;
   }

   int I3DBuilder::splittableTotal(const Dimensions::Transform::Id transId, Splitting::Locations::Id flag)
   {
      auto spTools = this->truncationTools(transId);
      auto n2D = this->dim(transId, Dimensions::Data::DAT2D);
      auto n3D = this->dim(transId, Dimensions::Data::DAT3D);

      auto nModes = spTools->totalModes(n2D, n3D, flag);

      // If none of the previous conditions were right
      if(nModes < 0)
      {
         throw std::logic_error("Tried to split in a unknown dimension for 3D spatial scheme builder");
      }

      return nModes;
   }

   void I3DBuilder::split(std::multimap<int,int>& modes, const std::vector<int>& id, const std::vector<int>& bins, const Dimensions::Transform::Id transId)
   {
      auto spTools = this->truncationTools(transId);
      std::vector<int> n0_, nN_;

      auto n2D = this->dim(transId, Dimensions::Data::DAT2D);
      auto n3D = this->dim(transId, Dimensions::Data::DAT3D);

      // Create list of modes
      spTools->buildBalancedMap(modes, n2D, n3D, id, bins, n0_, nN_);
   }

   void I3DBuilder::setDimensions()
   {
      //
      // Set transform space sizes
      //
      ArrayI traSize(3);
      traSize(0) = this->mesher().nSpec1D();
      traSize(1) = this->mesher().nSpec2D();
      traSize(2) = this->mesher().nSpec3D();
      this->setTransformSpace(traSize);

      //
      // Initialise spectral space
      //

      // Initialise forward dimension of first transform
      this->setDimension(this->mesher().nSpec1D(), Dimensions::Transform::SPECTRAL, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(this->mesher().nSpec1D(), Dimensions::Transform::SPECTRAL, Dimensions::Data::DATB1D);

      // Initialise second dimension of first transform
      this->setDimension(this->mesher().nSpec3D(), Dimensions::Transform::SPECTRAL, Dimensions::Data::DAT2D);

      // Initialise third dimension of first transform
      this->setDimension(this->mesher().nSpec2D(), Dimensions::Transform::SPECTRAL, Dimensions::Data::DAT3D);

      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->setDimension(this->mesher().nPhys1D(), Dimensions::Transform::TRA1D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(this->mesher().nDealias1D(), Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D);

      // Initialise second dimension of first transform
      this->setDimension(this->mesher().nSpec3D(), Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

      // Initialise third dimension of first transform
      this->setDimension(this->mesher().nSpec2D(), Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);

      //
      // Initialise second transform
      //

      // Initialise forward dimension of second transform
      this->setDimension(this->mesher().nPhys2D(), Dimensions::Transform::TRA2D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of second transform
      this->setDimension(this->mesher().nDealias2D(), Dimensions::Transform::TRA2D, Dimensions::Data::DATB1D);

      // Initialise second dimension of second transform
      this->setDimension(this->mesher().nPhys1D(), Dimensions::Transform::TRA2D, Dimensions::Data::DAT2D);

      // Initialise third dimension of second transform
      this->setDimension(this->mesher().nSpec3D(), Dimensions::Transform::TRA2D, Dimensions::Data::DAT3D);

      //
      // Initialise third transform
      //

      // Initialise forward dimension of third transform
      this->setDimension(this->mesher().nPhys3D(), Dimensions::Transform::TRA3D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of third transform
      this->setDimension(this->mesher().nDealias3D(), Dimensions::Transform::TRA3D, Dimensions::Data::DATB1D);

      // Initialise second dimension of third transform
      this->setDimension(this->mesher().nPhys2D(), Dimensions::Transform::TRA3D, Dimensions::Data::DAT2D);

      // Initialise third dimension of third transform
      this->setDimension(this->mesher().nPhys1D(), Dimensions::Transform::TRA3D, Dimensions::Data::DAT3D);
   }

   void I3DBuilder::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);

      // Set second transform cost
      this->setCost(1.0, Dimensions::Transform::TRA2D);

      // Set third transform cost
      this->setCost(1.0, Dimensions::Transform::TRA3D);
   }

   void I3DBuilder::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);

      // Set second transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA2D);

      // Set third transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA3D);
   }

   void I3DBuilder::setMemoryScore()
   {
      // Set first transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA1D);

      // Set second transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA2D);

      // Set third transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA3D);
   }

} // SpatialScheme
} // QuICC
