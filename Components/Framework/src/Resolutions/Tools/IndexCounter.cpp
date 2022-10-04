/** 
 * @file IndexCounter.cpp
 * @brief Source of base class for the index counters
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"

// Project includes
//

namespace QuICC {

   IndexCounter::IndexCounter()
   {
   }

   IndexCounter::~IndexCounter()
   {
   }

   std::vector<int> IndexCounter::makeVKey(const Dimensions::Transform::Id id, const int i, const int j, const int k) const
   {
      std::vector<int> key;

      if(id == Dimensions::Transform::TRA1D || id == Dimensions::Transform::SPECTRAL)
      {
         key.push_back(i);
         key.push_back(k);
         key.push_back(j);
      }
      else if(id == Dimensions::Transform::TRA2D)
      {
         key.push_back(j);
         key.push_back(i);
         key.push_back(k);
      }
      else if(id == Dimensions::Transform::TRA3D)
      {
         key.push_back(k);
         key.push_back(j);
         key.push_back(i);
      }
      else
      {
         throw std::logic_error("Unknown ID used to generate VKey");
      }

      return key;
   }

   std::vector<int> IndexCounter::makeVKey(const Dimensions::Transform::Id id, const int i, const int j) const
   {
      std::vector<int> key;

      if(id == Dimensions::Transform::TRA1D || id == Dimensions::Transform::SPECTRAL)
      {
         key.push_back(i);
         key.push_back(j);
      }
      else if(id == Dimensions::Transform::TRA2D)
      {
         key.push_back(j);
         key.push_back(i);
      }
      else if(id == Dimensions::Transform::TRA3D)
      {
         throw std::logic_error("Tried to use 2D keymaker on third dimension");
      }
      else
      {
         throw std::logic_error("Unknown ID used to generate VKey");
      }

      return key;
   }

   std::tuple<int,int,int> IndexCounter::makeKey(const Dimensions::Transform::Id id, const int i, const int j, const int k) const
   {
      std::tuple<int,int,int> key;

      if(id == Dimensions::Transform::TRA1D || id == Dimensions::Transform::SPECTRAL)
      {
         key = std::make_tuple(i, k, j);
      }
      else if(id == Dimensions::Transform::TRA2D)
      {
         key = std::make_tuple(j, i, k);
      }
      else if(id == Dimensions::Transform::TRA3D)
      {
         key = std::make_tuple(k, j, i);
      }
      else
      {
         throw std::logic_error("Unknown ID used to generate Key");
      }

      return key;
   }

   std::pair<int,int> IndexCounter::makeKey(const Dimensions::Transform::Id id, const int i, const int j) const
   {
      std::pair<int,int> key;

      if(id == Dimensions::Transform::TRA1D || id == Dimensions::Transform::SPECTRAL)
      {
         key = std::make_pair(i, j);
      }
      else if(id == Dimensions::Transform::TRA2D)
      {
         key = std::make_pair(j, i);
      }
      else if(id == Dimensions::Transform::TRA3D)
      {
         throw std::logic_error("Tried to use 2D keymaker on third dimension");
      }
      else
      {
         throw std::logic_error("Unknown ID used to generate Key");
      }

      return key;
   }

}
