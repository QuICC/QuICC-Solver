/** 
 * @file IRegular3DBuilder.cpp
 * @brief Source of the regular 3D scheme builder
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/3D/IRegular3DBuilder.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/Tools/Regular.hpp"

namespace QuICC {

namespace SpatialScheme {

   IRegular3DBuilder::IRegular3DBuilder(const ArrayI& dim, const GridPurpose::Id purpose)
      : IBuilder(dim.size(), purpose), mI(dim(0)), mJ(dim(1)), mK(dim(2))
   {
      assert(dim.size() == 3);
   }

   IRegular3DBuilder::~IRegular3DBuilder()
   {
   }

   ArrayI IRegular3DBuilder::resolution() const
   {
      ArrayI space(this->dims());
      space << this->mI, this->mJ, this->mK;

      return space;
   }

   int IRegular3DBuilder::fillIndexes(const Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, const Splitting::Locations::Id flag)
   {
      // Safety assertions for default values
      assert( (id.size() == 0) || (bins.size() > 0) );
      assert( id.size() == bins.size() );
      assert( n0.size() == nN.size() );
      assert( (bins.size() == 0) || (flag != Splitting::Locations::NONE) );

      // Make sure we start with empty indexes
      fwd1D.clear();
      bwd1D.clear();
      idx2D.clear();

      // Multimap for the modes
      std::multimap<int,int> modes;

      // No splitting
      if(flag == Splitting::Locations::NONE)
      {
         this->splitSerial(modes, transId);

      // Splitting is on first transform
      } else if(flag == Splitting::Locations::FIRST)
      {
         this->splitSingle1D(modes, n0, nN, transId);

      // Splitting is on second transform
      } else if(flag == Splitting::Locations::SECOND)
      {
         this->splitSingle2D(modes, n0, nN, transId);

      // Splitting is on both transforms
      } else if(flag == Splitting::Locations::BOTH)
      {
         this->splitTubular(modes, n0, nN, transId);

      // Splitting is on slowest index on first transforms
      } else if(flag == Splitting::Locations::COUPLED2D)
      {
         this->splitCoupled2D(modes, n0, nN, transId);
      }

      // Fill indexes for 2D and 3D
      Tools::Regular::fillIndexes2D3D(idx2D, idx3D, modes);

      // Fill indexes for 1D
      Tools::Regular::fillIndexes1D(fwd1D, bwd1D, idx3D, this->dim(transId, Dimensions::Data::DATF1D), this->dim(transId, Dimensions::Data::DATB1D));

      // Set status (0 for success, 1 for failure)
      int status = 0;
      if(modes.size() == 0)
      {
         status = 1;
      }

      return status;
   }

   int IRegular3DBuilder::splittableTotal(const Dimensions::Transform::Id transId, Splitting::Locations::Id flag)
   {
      // Splittable size for first transform splitting
      if(flag == Splitting::Locations::FIRST)
      {
         // Get total size for first transform
         if(transId == Dimensions::Transform::TRA1D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

         // Get total size for second transform
         } else if(transId == Dimensions::Transform::TRA2D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D);

         // Get total size for third transform
         } else if(transId == Dimensions::Transform::TRA3D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);
         }

      // Splittable size for second transform splitting
      } else if(flag == Splitting::Locations::SECOND)
      {
         // Get total size for first transform
         if(transId == Dimensions::Transform::TRA1D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D);
         // Get total size for second transform
         } else if(transId == Dimensions::Transform::TRA2D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);

         // Get total size for third transform
         } else if(transId == Dimensions::Transform::TRA3D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);
         }

      // Splittable size for both transforms splitting
      } else if(flag == Splitting::Locations::BOTH)
      {
         // Get total size for first transform
         if(transId == Dimensions::Transform::TRA1D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);

         // Get total size for second transform
         } else if(transId == Dimensions::Transform::TRA2D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);

         // Get total size for third transform
         } else if(transId == Dimensions::Transform::TRA3D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D);
         }

      // Splittable size for coupled 2d matrices first transforms splitting
      } else if(flag == Splitting::Locations::COUPLED2D)
      {
         // Get total size for first transform
         if(transId == Dimensions::Transform::TRA1D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);

         // Get total size for second transform
         } else if(transId == Dimensions::Transform::TRA2D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D);

         // Get total size for third transform
         } else if(transId == Dimensions::Transform::TRA3D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);
         }
      }
      
      // If none of the previous conditions were right
      throw std::logic_error("Tried to split in a unknown dimension for 3D regular case");

      return -1;
   }

   void IRegular3DBuilder::splitSerial(std::multimap<int,int>& modes, const Dimensions::Transform::Id transId)
   {
      int k0 = 0;
      int kN = this->dim(transId, Dimensions::Data::DAT3D);
      ArrayI j0 = ArrayI::Zero(kN);
      ArrayI jN = ArrayI::Constant(kN, this->dim(transId, Dimensions::Data::DAT2D));
      int c0 = 0;
      int cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

      // Generate map for regular indexes
      Tools::Regular::buildMap(modes, k0, kN, j0, jN, c0, cN);
   }

   void IRegular3DBuilder::splitSingle1D(std::multimap<int,int>& modes, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId)
   {
      int k0 = -1;
      int kN = -1;
      ArrayI j0, jN;
      int c0 = -1;
      int cN = -1;

      // Create index list for first transform
      if(transId == Dimensions::Transform::TRA1D)
      {
         k0 = 0;
         kN = this->dim(transId, Dimensions::Data::DAT3D);
         j0 = ArrayI::Zero(kN);
         jN = ArrayI::Constant(kN, this->dim(transId, Dimensions::Data::DAT2D));
         c0 = n0(0);
         cN = n0(0) + nN(0);

      // Create index list for second transform
      } else if(transId == Dimensions::Transform::TRA2D)
      {
         k0 = 0;
         kN = this->dim(transId, Dimensions::Data::DAT3D);
         j0 = ArrayI::Constant(kN, n0(0));
         jN = ArrayI::Constant(kN, nN(0));
         c0 = 0;
         cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

      // Create index list for third transform
      } else if(transId == Dimensions::Transform::TRA3D)
      {
         k0 = n0(0);
         kN = nN(0);
         j0 = ArrayI::Zero(kN);
         jN = ArrayI::Constant(kN, this->dim(transId, Dimensions::Data::DAT2D));
         c0 = 0;
         cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);
      }

      // Generate map for regular indexes
      Tools::Regular::buildMap(modes, k0, kN, j0, jN, c0, cN);
   }

   void IRegular3DBuilder::splitSingle2D(std::multimap<int,int>& modes, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId)
   {
      int k0 = -1;
      int kN = -1;
      ArrayI j0, jN;
      int c0 = -1;
      int cN = -1;

      // Create index list for first transform
      if(transId == Dimensions::Transform::TRA1D)
      {
         k0 = 0;
         kN = this->dim(transId, Dimensions::Data::DAT3D);
         j0 = ArrayI::Constant(kN, n0(0));
         jN = ArrayI::Constant(kN, nN(0));
         c0 = 0;
         cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

      // Create index list for second transform
      } else if(transId == Dimensions::Transform::TRA2D)
      {
         k0 = n0(0);
         kN = nN(0);
         j0 = ArrayI::Zero(kN);
         jN = ArrayI::Constant(kN, this->dim(transId, Dimensions::Data::DAT2D));
         c0 = 0;
         cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

      // Create index list for third transform
      } else if(transId == Dimensions::Transform::TRA3D)
      {
         k0 = 0;
         kN = this->dim(transId, Dimensions::Data::DAT3D);
         j0 = ArrayI::Zero(kN);
         jN = ArrayI::Constant(kN, this->dim(transId, Dimensions::Data::DAT2D));
         c0 = n0(0);
         cN = n0(0) + nN(0);
      }

      // Generate map for regular indexes
      Tools::Regular::buildMap(modes, k0, kN, j0, jN, c0, cN);
   }

   void IRegular3DBuilder::splitCoupled2D(std::multimap<int,int>& modes, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId)
   {
      int k0 = -1;
      int kN = -1;
      ArrayI j0, jN;
      int c0 = -1;
      int cN = -1;

      // Create index list for first transform
      if(transId == Dimensions::Transform::TRA1D)
      {
            k0 = n0(0);
            kN = nN(0);
            j0 = ArrayI::Zero(kN);
            jN = ArrayI::Constant(kN, this->dim(transId, Dimensions::Data::DAT2D));
            c0 = 0;
            cN = this->dim(transId, Dimensions::Data::DAT2D)*kN;

      // Create index list for second transform
      } else if(transId == Dimensions::Transform::TRA2D)
      {
            k0 = 0;
            kN = this->dim(transId, Dimensions::Data::DAT3D);
            j0 = ArrayI::Constant(kN, n0(0));
            jN = ArrayI::Constant(kN, nN(0));
            c0 = 0;
            cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

      // Create index list for third transform
      } else if(transId == Dimensions::Transform::TRA3D)
      {
            k0 = n0(0);
            kN = nN(0);
            j0 = ArrayI::Zero(kN);
            jN = ArrayI::Constant(kN, this->dim(transId, Dimensions::Data::DAT2D));
            c0 = 0;
            cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);
      }

      // Generate map for regular indexes
      Tools::Regular::buildMap(modes, k0, kN, j0, jN, c0, cN);
   }

   void IRegular3DBuilder::splitTubular(std::multimap<int,int>& modes, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId)
   {
      int k0 = -1;
      int kN = -1;
      ArrayI j0, jN;
      int c0 = -1;
      int cN = -1;

      // Create index list for first transform
      if(transId == Dimensions::Transform::TRA1D)
      {
         k0 = n0(0);
         kN = nN(0);
         j0 = n0.tail(kN);
         jN = nN.tail(kN);
         c0 = 0;
         cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

      // Create index list for second transform
      } else if(transId == Dimensions::Transform::TRA2D)
      {
         k0 = n0(0);
         kN = nN(0);
         j0 = ArrayI::Constant(kN, n0(1));
         jN = ArrayI::Constant(kN, nN(1));
         c0 = 0;
         cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

      // Create index list for third transform
      } else if(transId == Dimensions::Transform::TRA3D)
      {
         k0 = n0(0);
         kN = nN(0);
         j0 = n0.tail(kN);
         jN = nN.tail(kN);
         c0 = 0;
         cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);
      }

      // Generate map for regular indexes
      Tools::Regular::buildMap(modes, k0, kN, j0, jN, c0, cN);
   }
}
}
