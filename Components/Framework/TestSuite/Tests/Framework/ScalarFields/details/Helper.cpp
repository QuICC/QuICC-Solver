#include <algorithm>
#include "Helper.hpp"

namespace details {

   std::shared_ptr<QuICC::Datatypes::ScalarFieldSetup> createSetup(const SetupType type, const std::size_t dim1D, const std::size_t dim3D, const std::vector<std::size_t>& idx3D)
   {
      auto spMeta = std::make_shared<QuICC::CscMetadata>();
      auto& meta = *spMeta;

      meta.global3D = dim3D;

      meta.ptr2D.reserve(dim3D + 1);
      meta.ptr2D.push_back(0);
      auto it3D = idx3D.cbegin();
      if(type == SetupType::UniformUp)
      {
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
         meta.global2D = std::max(meta.global2D, meta.ptr2D.at(i+1) - meta.ptr2D.at(i));
      }

      for(std::uint32_t i = 0; i < meta.ptr2D.size()-1; i++)
      {
         auto sze = meta.ptr2D.at(i + 1) - meta.ptr2D.at(i);
         for(std::uint32_t j = 0; j < sze; j++)
         {
            meta.idx2D.push_back(j);
         }
      }

      auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();

      auto spSetup = std::make_shared<QuICC::Datatypes::ScalarFieldSetup>(spMeta, mem);

      return spSetup;
   }
}

