/**
 * @file FieldProbeWriter.cpp
 * @brief Source of the implementation of the physical space field probe
 */

// System includes
//
#include <iomanip>
#include <limits>
#include <locale>
#include <stdexcept>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Io/Variable/FieldProbeWriter.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "Types/Math.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Io/Variable/Tags/FieldProbe.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   FieldProbeWriter::FieldProbeWriter(const std::string& prefix, const std::string& type, const std::vector<MHDFloat>& position)
      : IVariableAsciiWriter(prefix + Tags::FieldProbe::BASENAME, Tags::FieldProbe::EXTENSION, prefix + Tags::FieldProbe::HEADER, type, Tags::FieldProbe::VERSION, Dimensions::Space::PHYSICAL, EXTEND), mPosition(position)
   {
   }

   void FieldProbeWriter::init()
   {
      if(this->mMesh.size() != this->mPosition.size())
      {
         throw std::logic_error("Mesh and position do not have the same dimensionality");
      }

      // Search for nearest grid point
      std::vector<MHDFloat> nearest;
      std::vector<int> indexes;
      for(std::size_t i = 0; i < this->mPosition.size(); i++)
      {
         const auto& coord = this->mPosition.at(i);
         const auto& grid = this->mMesh.at(i);
         auto dist = std::numeric_limits<MHDFloat>::max();
         int idx = -1;
         for(int j = 0; j < grid.size(); j++)
         {
            auto g = grid(j);
            auto d = std::abs(coord - g);
            if(d < dist)
            {
               dist = d;
               idx = j;
            }
         }
         if(idx == -1)
         {
            throw std::logic_error("Coordinate was out of range");
         }

         nearest.push_back(grid(idx));
         indexes.push_back(idx);
      }

      // Update position with nearest grid point
      this->mPosition = nearest;

      //Check if first coordinate is stored on local rank
      const auto& tRes = *this->res().cpu()->dim(Dimensions::Transform::TRA3D);
      bool found = false;
      const auto& coordIdx = indexes.at(0);
      for(auto k = 0;k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
      {
         auto k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
         if(k_ == coordIdx)
         {
            found = true;
            this->mIndexes.push_back(k);
            break;
         }
         else if(k_ > coordIdx)
         {
            break;
         }
      }

      //Check if second coordinate is stored on local rank
      if(found)
      {
         found = false;
         const auto k = this->mIndexes.front();
         const auto& coordIdx = indexes.at(1);
         for(auto j = 0;j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
         {
            auto j_ = tRes.idx<Dimensions::Data::DAT2D>(j,k);
            if(j_ == coordIdx)
            {
               found = true;
               this->mIndexes.push_back(j);
               break;
            }
            else if(j_ > coordIdx)
            {
               break;
            }
         }
      }

      //Check if third coordinate is stored on local rank
      if(found)
      {
         found = false;
         const auto k = this->mIndexes.front();
         const auto j = this->mIndexes.at(1);
         const auto& coordIdx = indexes.at(2);
         for(auto i = 0;i < tRes.dim<Dimensions::Data::DATF1D>(j, k); i++)
         {
            auto i_ = tRes.idx<Dimensions::Data::DATF1D>(i,j,k);
            if(i_ == coordIdx)
            {
               found = true;
               this->mIndexes.push_back(i);
               break;
            }
            else if(i_ > coordIdx)
            {
               break;
            }
         }
      }

      // If not present, clear position and indexes
      if(!found)
      {
         this->mIndexes.clear();
      }

      IVariableAsciiWriter::init();
   }

   void FieldProbeWriter::writeContent()
   {
      auto sRange = this->scalarRange();
      auto vRange = this->vectorRange();
      assert(std::distance(sRange.first, sRange.second) + std::distance(vRange.first, vRange.second) == 1);

      auto isScalar = (std::distance(sRange.first, sRange.second) == 1);
      auto isVector = (std::distance(vRange.first, vRange.second) == 1);

      if(this->mIndexes.size() > 0)
      {
         this->mValue.clear();

         auto i = this->mIndexes.back();
         auto j = this->mIndexes.at(1);
         auto k = this->mIndexes.front();

         if(isScalar)
         {
            auto v = std::visit([&](auto&& p){return p->dom(0).phys().comp(FieldComponents::Physical::SCALAR).point(i,j,k);}, sRange.first->second);
            this->mValue.push_back(v);
         }

         if(isVector)
         {
            auto v = std::visit([&](auto&& p){return p->dom(0).phys().comp(this->res().sim().ss().physical().ONE()).point(i,j,k);}, vRange.first->second);
            this->mValue.push_back(v);
            v = std::visit([&](auto&& p){return p->dom(0).phys().comp(this->res().sim().ss().physical().TWO()).point(i,j,k);}, vRange.first->second);
            this->mValue.push_back(v);
            v = std::visit([&](auto&& p){return p->dom(0).phys().comp(this->res().sim().ss().physical().THREE()).point(i,j,k);}, vRange.first->second);
            this->mValue.push_back(v);
         }
      } else
      {
         if(isScalar)
         {
            this->mValue = std::vector<MHDFloat>(1, 0);
         }

         if(isVector)
         {
            this->mValue = std::vector<MHDFloat>(3, 0);
         }
      }

      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mValue.data(), this->mValue.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      using Tools::Formatter::ioFW;
      int ioPrec = 14;

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::scientific;
         this->mFile << std::setprecision(ioPrec) << ioFW(ioPrec) << this->mTime << "\t" << ioFW(ioPrec);
         for(auto coord: this->mPosition)
         {
            this->mFile << "\t" << coord;
         }
         for(auto v: this->mValue)
         {
            this->mFile << "\t" << v;
         }
         this->mFile << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      for(auto v: this->mValue)
      {
         if(std::isnan(v))
         {
            QuICCEnv().abort("Field proble value is NaN!");
         }
      }
   }

} // Variable
} // Io
} // QuICC
