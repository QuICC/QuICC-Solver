#include <algorithm>
#include <random>
#include <cstdint>
#include <set>
#include "Helper.hpp"

namespace details {

   std::shared_ptr<QuICC::Datatypes::ScalarFieldSetup> createSetup(std::size_t& dim1D, std::size_t& dim2D, std::size_t& dim3D, std::vector<std::vector<std::size_t>>& idx2D, std::vector<std::size_t>& idx3D, const SetupType type, const std::uint32_t variant)
   {
      std::vector<int> idxList3D;
      std::vector<int> idxList2D;
      std::set<int> sorter;
      std::random_device rd;
      std::mt19937 g(rd());

      if(variant < 10)
      {
         dim1D = 16;
         dim2D = 16;
         dim3D = 16;
      }
      else if(variant < 20)
      {
         dim1D = 30;
         dim2D = 30;
         dim3D = 30;
      }
      else if(variant < 30)
      {
         dim1D = 64;
         dim2D = 64;
         dim3D = 64;
      }
      else if(variant < 40)
      {
         dim1D = 100;
         dim2D = 100;
         dim3D = 100;
      }
      else if(variant < 50)
      {
         dim1D = 256;
         dim2D = 256;
         dim3D = 256;
      }
      else
      {
         throw std::logic_error("Unknown variant");
      }

      idxList2D.resize(dim2D);
      std::iota(idxList2D.begin(), idxList2D.end(), 0);
      std::shuffle(idxList2D.begin(), idxList2D.end(), g);

      idxList3D.resize(dim3D);
      std::iota(idxList3D.begin(), idxList3D.end(), 0);
      std::shuffle(idxList3D.begin(), idxList3D.end(), g);
      std::uniform_int_distribution<> n2D(3, static_cast<std::uint32_t>(0.5*dim2D));
      std::uniform_int_distribution<> n3D(3, static_cast<std::uint32_t>(0.5*dim3D));

      // Get 3D indexes
      idx3D.clear();
      std::uint32_t sze3D = n3D(g);
      std::copy_n(idxList3D.begin(), sze3D, std::inserter(sorter, sorter.end()));
      assert(sorter.size() == sze3D);
      idx3D.reserve(sze3D);
      std::copy_n(sorter.begin(), sze3D, std::back_inserter(idx3D));
      assert(idx3D.size() == sze3D);
      auto it3D = idx3D.cbegin();

      auto spMeta = std::make_shared<QuICC::CscMetadata>();
      auto& meta = *spMeta;

      // Set global indexes
      meta.global1D = dim1D;
      meta.global2D = dim2D;
      meta.global3D = dim3D;

      meta.ptr2D.reserve(dim3D + 1);
      meta.ptr2D.push_back(0);
      if(type == SetupType::UniformUp)
      {
         for(std::size_t i = 0; i < dim3D; i++)
         {
            if(it3D != idx3D.end() && i == *it3D)
            {
               std::uint32_t sze2D = n2D(g);
               meta.ptr2D.push_back(meta.ptr2D.back() + sze2D);
               for(std::size_t j = 0; j < sze2D; j++)
               {
                  meta.dim1D.push_back(dim1D);
               }
               std::shuffle(idxList2D.begin(), idxList2D.end(), g);
               sorter.clear();
               std::copy_n(idxList2D.begin(), sze2D, std::inserter(sorter, sorter.end()));
               assert(sorter.size() == sze2D);
               std::copy_n(sorter.begin(), sze2D, std::back_inserter(meta.idx2D));
               idx2D.push_back({});
               idx2D.back().reserve(sze2D);
               std::copy_n(sorter.begin(), sze2D, std::back_inserter(idx2D.back()));
               assert(idx2D.back().size() == sze2D);
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

      auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();

      auto spSetup = std::make_shared<QuICC::Datatypes::ScalarFieldSetup>(spMeta, mem);

      return spSetup;
   }
}

