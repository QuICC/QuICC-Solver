/** 
 * @file IRegular1DBuilder.cpp
 * @brief Source of the generic regular 1D scheme implementation
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/1D/IRegular1DBuilder.hpp"

// Project includes
//

namespace QuICC {

namespace SpatialScheme {

   IRegular1DBuilder::IRegular1DBuilder(const ArrayI& dim, const GridPurpose::Id purpose)
      : IBuilder(dim.size(), purpose), mI(dim(0))
   {
      assert(dim.size() == 1);
   }

   IRegular1DBuilder::~IRegular1DBuilder()
   {
   }

   ArrayI IRegular1DBuilder::resolution() const
   {
      ArrayI space(this->dims());
      space << this->mI;

      return space;
   }

   int IRegular1DBuilder::fillIndexes(Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, const Splitting::Locations::Id flag)
   {
      // Assert for right transform (1D case)
      assert(transId == Dimensions::Transform::TRA1D);

      // Safety assertions for default values
      assert( bins.size() == n0.size() );
      assert( n0.size() == nN.size() );
      assert( (bins.size() == 0) || (flag != Splitting::Locations::NONE) );

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

}
}
