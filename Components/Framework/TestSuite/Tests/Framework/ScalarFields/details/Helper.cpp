#include <algorithm>
#include <cstdint>
#include "Helper.hpp"

namespace details {

   std::shared_ptr<QuICC::Datatypes::ScalarFieldSetup> createSetup(std::size_t& dim1D, std::size_t& dim3D, std::vector<std::size_t>& idx3D, const SetupType type, const std::uint32_t variant)
   {
      if(variant == 0)
      {
         dim3D = 5;
         idx3D = {0, 1, 2, 3, 4};
      }
      else if(variant == 1)
      {
         dim3D = 10;
         idx3D = {2, 3, 5, 6, 8};
      }
      else
      {
         throw std::logic_error("Unknown variant");
      }

      auto spMeta = std::make_shared<QuICC::CscMetadata>();
      auto& meta = *spMeta;

      meta.global3D = dim3D;

      meta.ptr2D.reserve(dim3D + 1);
      meta.ptr2D.push_back(0);
      auto it3D = idx3D.cbegin();
      if(type == SetupType::UniformUp)
      {
         dim1D = 2*dim3D;
         for(std::size_t i = 0; i < dim3D; i++)
         {
            if(i == *it3D)
            {
               std::uint32_t sze = (i + 2);
               meta.ptr2D.push_back(meta.ptr2D.back() + sze);
               for(std::size_t j = 0; j < sze; j++)
               {
                  meta.dim1D.push_back(dim1D);
               }
               it3D++;
            }
            else
            {
               meta.ptr2D.push_back(meta.ptr2D.back());
            }
         }
      }
      else if(type == SetupType::UniformDown)
      {
         dim1D = 2*dim3D;
         for(std::size_t i = 0; i < dim3D; i++)
         {
            if(i == *it3D)
            {
               std::uint32_t sze = (dim1D - i);
               meta.ptr2D.push_back(meta.ptr2D.back() + sze);
               for(std::size_t j = 0; j < sze; j++)
               {
                  meta.dim1D.push_back(dim1D);
               }
               it3D++;
            }
            else
            {
               meta.ptr2D.push_back(meta.ptr2D.back());
            }
         }
      }
      else if(type == SetupType::TriangularUp)
      {
         dim1D = 2;
         for(std::size_t i = 0; i < dim3D; i++)
         {
            if(i == *it3D)
            {
               std::uint32_t sze = (i + 1);
               meta.ptr2D.push_back(meta.ptr2D.back() + sze);
               for(std::size_t j = 0; j < sze; j++)
               {
                  meta.dim1D.push_back(i + dim1D);
               }
               it3D++;
            }
            else
            {
               meta.ptr2D.push_back(meta.ptr2D.back());
            }
         }
      }
      else if(type == SetupType::TriangularDown)
      {
         dim1D = 3*dim3D;
         for(std::size_t i = 0; i < dim3D; i++)
         {
            if(i == *it3D)
            {
               std::uint32_t sze = (dim1D - 2*i);
               meta.ptr2D.push_back(meta.ptr2D.back() + sze);
               for(std::size_t j = 0; j < sze; j++)
               {
                  meta.dim1D.push_back(dim1D - i);
               }
               it3D++;
            }
            else
            {
               meta.ptr2D.push_back(meta.ptr2D.back());
            }
         }
      }
      else
      {
         throw std::logic_error("Unknown setup type");
      }

      // Global dimension
      meta.global1D = *std::max_element(meta.dim1D.begin(), meta.dim1D.end());
      meta.global2D = 0;
      for(std::size_t i = 0; i < meta.ptr2D.size() - 1; i++)
      {
         meta.global2D = std::max(meta.global2D, meta.ptr2D.at(i+1) - meta.ptr2D.at(i) + 1);
      }

      for(std::uint32_t i = 0; i < meta.ptr2D.size()-1; i++)
      {
         auto sze = meta.ptr2D.at(i + 1) - meta.ptr2D.at(i);
         for(std::uint32_t j = 0; j < sze; j++)
         {
            // Shift by one to trigger error if local indexes are wrong
            meta.idx2D.push_back(j+1);
         }
      }

      auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();

      auto spSetup = std::make_shared<QuICC::Datatypes::ScalarFieldSetup>(spMeta, mem);

      return spSetup;
   }
}

