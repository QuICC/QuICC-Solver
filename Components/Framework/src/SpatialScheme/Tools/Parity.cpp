/** 
 * @file Parity.cpp
 * @brief Source of the tools for parity splitting
 */

// System includes
//
#include <set>
#include <map>
#include <tuple>

// Project includes
//
#include "QuICC/SpatialScheme/Tools/Parity.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   void Parity::splitParityL(SharedResolution spRes, const Dimensions::Transform::Id traId, ArrayI& blockSize, MatrixI& evenBlocks, MatrixI& oddBlocks)
   {
      // Get number of transforms
      std::vector<std::tuple<int,int,int> > even;
      std::vector<std::tuple<int,int,int> > odd;
      int idx = 0;
      int previous = -1;
      for(int i = 0; i < spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         if(spRes->cpu()->dim(traId)->idx<Dimensions::Data::DAT3D>(i)%2 == 0)
         {
            if(previous == 0)
            {
               std::get<1>(even.back()) += spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i);
            } else
            {
               even.push_back(std::make_tuple(idx, spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i), spRes->cpu()->dim(traId)->idx<Dimensions::Data::DAT3D>(i)));
            }

            idx += spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i);

            previous = 0;
         } else
         {
            if(previous == 1)
            {
               std::get<1>(odd.back()) += spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i);
            } else
            {
               odd.push_back(std::make_tuple(idx, spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i),spRes->cpu()->dim(traId)->idx<Dimensions::Data::DAT3D>(i)));
            }

            idx += spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i);

            previous = 1;
         }
      }

      evenBlocks.resize(even.size(), 3);
      oddBlocks.resize(odd.size(), 3);

      int i = 0;
      for(auto it = even.cbegin(); it != even.cend(); ++it, ++i)
      {
         evenBlocks(i,0) = std::get<0>(*it);
         evenBlocks(i,1) = std::get<1>(*it);
         evenBlocks(i,2) = std::get<2>(*it);
      }

      i = 0;
      for(auto it = odd.cbegin(); it != odd.cend(); ++it, ++i)
      {
         oddBlocks(i,0) = std::get<0>(*it);
         oddBlocks(i,1) = std::get<1>(*it);
         oddBlocks(i,2) = std::get<2>(*it);
      }

      blockSize.resize(2);
      blockSize(0) = evenBlocks.col(1).sum();
      blockSize(1) = oddBlocks.col(1).sum();
   }

   void Parity::splitParityM(SharedResolution spRes, const Dimensions::Transform::Id traId, ArrayI& blockSize, MatrixI& evenBlocks, MatrixI& oddBlocks)
   {
      // Get number of transforms
      std::vector<std::tuple<int,int,int> > even;
      std::vector<std::tuple<int,int,int> > odd;
      int idx = 0;
      int previous = -1;
      for(int i = 0; i < spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         for(int j = 0; j < spRes->cpu()->dim(traId)->dim<Dimensions::Data::DAT2D>(i); j++)
         {
            if(spRes->cpu()->dim(traId)->idx<Dimensions::Data::DAT2D>(j,i)%2 == 0)
            {
               if(previous == 0)
               {
                  std::get<1>(even.back()) += 1;
               } else
               {
                  even.push_back(std::make_tuple(idx, 1, spRes->cpu()->dim(traId)->idx<Dimensions::Data::DAT2D>(j,i)));
               }

               idx += 1;

               previous = 0;
            }
            else
            {
               if(previous == 1)
               {
                  std::get<1>(odd.back()) += 1;
               } else
               {
                  odd.push_back(std::make_tuple(idx, 1, spRes->cpu()->dim(traId)->idx<Dimensions::Data::DAT2D>(j,i)));
               }

               idx += 1;

               previous = 1;
            }
         }
      }

      evenBlocks.resize(even.size(), 3);
      oddBlocks.resize(odd.size(), 3);

      int i = 0;
      for(auto it = even.cbegin(); it != even.cend(); ++it, ++i)
      {
         evenBlocks(i,0) = std::get<0>(*it);
         evenBlocks(i,1) = std::get<1>(*it);
         evenBlocks(i,2) = std::get<2>(*it);
      }

      i = 0;
      for(auto it = odd.cbegin(); it != odd.cend(); ++it, ++i)
      {
         oddBlocks(i,0) = std::get<0>(*it);
         oddBlocks(i,1) = std::get<1>(*it);
         oddBlocks(i,2) = std::get<2>(*it);
      }

      blockSize.resize(2);
      blockSize(0) = evenBlocks.col(1).sum();
      blockSize(1) = oddBlocks.col(1).sum();
   }

} // Tools
} // SpatialScheme
} // QuICC
