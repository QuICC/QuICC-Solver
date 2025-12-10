#include <algorithm>
#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

#include "details/Helper.hpp"

TEST_CASE("Uniform increasing ScalarField setup test", "[Datatypes::ScalarFieldSetup_uniformUp]")
{
   std::size_t dim3D;
   std::size_t dim1D;
   std::vector<std::size_t> idx3D;
   auto spSetup = details::createSetup(dim1D, dim3D, idx3D, details::SetupType::UniformUp);

   auto dfct = [](const std::size_t k){return k + 2;};
   auto sze = [&](const std::size_t& a, const auto& b){return a + dfct(b);};
   auto acc = [&](const std::size_t n){return static_cast<int>(std::accumulate(idx3D.begin(), idx3D.begin()+n, 0, sze));};

   int d3D = static_cast<int>(idx3D.size());
   int d1D = static_cast<int>(dim1D);

   CHECK( spSetup->nBlock() ==  d3D);
   CHECK( spSetup->dataRows() == d1D );
   CHECK( spSetup->dataCols() == acc(idx3D.size()) );
   CHECK( spSetup->blockIdx(1) == acc(1) );
   CHECK( spSetup->blockRows(2) == d1D );
   CHECK( spSetup->blockCols(3) == static_cast<int>(dfct(idx3D.at(3))) );
   CHECK( spSetup->colIdx(1) == 1 );
   CHECK( spSetup->colIdx(1, 0) == 1 );
   CHECK( spSetup->colIdx(2, 1) == acc(1) + 2 );
   CHECK( spSetup->colIdx(2, 3) == acc(3) + 2 );
}

TEST_CASE("Uniform increasing ScalarField setup global test", "[Datatypes::ScalarFieldSetup_global_uniformUp]")
{
   std::size_t dim3D;
   std::size_t dim1D;
   std::vector<std::size_t> idx3D;
   auto spSetup = details::createSetup(dim1D, dim3D, idx3D, details::SetupType::UniformUp, 1);

   auto dfct = [](const std::size_t k){return k + 2;};
   auto sze = [&](const std::size_t& a, const auto& b){return a + dfct(b);};
   auto acc = [&](const std::size_t n){return std::accumulate(idx3D.begin(), idx3D.begin()+n, 0, sze);};

   int d3D = static_cast<int>(idx3D.size());
   int d1D = static_cast<int>(dim1D);

   CHECK( spSetup->nBlock() == d3D );
   CHECK( spSetup->dataRows() == d1D );
   CHECK( spSetup->dataCols() == acc(idx3D.size()) );
   CHECK( spSetup->blockIdx(1) == acc(1) );
   CHECK( spSetup->blockRows(2) == d1D );
   CHECK( spSetup->blockCols(3) == static_cast<int>(dfct(idx3D.at(3))) );
   CHECK( spSetup->colIdx(1) == 1 );
   CHECK( spSetup->colIdx(1, 0) == 1 );
   CHECK( spSetup->colIdx(2, 1) == acc(1) + 2 );
   CHECK( spSetup->colIdx(2, 3) == acc(3) + 2 );
}

TEST_CASE("Uniform increasing ScalarField setup view test", "[Datatypes::ScalarFieldSetup_view_uniformUp]")
{
   std::size_t dim3D;
   std::size_t dim1D;
   std::vector<std::size_t> idx3D;
   auto spSetup = details::createSetup(dim1D, dim3D, idx3D, details::SetupType::UniformUp);
   auto&& meta = *spSetup->viewMeta();

   auto dfct = [](const std::size_t k){return k + 2;};
   auto sze = [&](const std::size_t& a, const auto& b){return a + dfct(b);};
   auto acc = [&](const std::size_t n){return std::accumulate(idx3D.begin(), idx3D.begin()+n, std::size_t(0), sze);};

   CHECK( meta.ptr2D.size() - 1 == idx3D.size() );
   CHECK( *std::max_element(meta.dim1D.begin(), meta.dim1D.end()) == dim1D );
   CHECK( *meta.ptr2D.rbegin() == acc(idx3D.size()) );
   CHECK( meta.ptr2D.at(1) == acc(1) );
   CHECK( meta.dim1D.at(2) == dim1D );
   CHECK( meta.ptr2D.at(4) - meta.ptr2D.at(3) == dfct(idx3D.at(3)) );
   CHECK( (1 + meta.ptr2D.at(0)) == 1 );
   CHECK( (1 + meta.ptr2D.at(0)) == 1 );
   CHECK( (2 + meta.ptr2D.at(1)) == acc(1) + 2 );
   CHECK( (2 + meta.ptr2D.at(3)) == acc(3) + 2 );
}

TEST_CASE("Uniform increasing ScalarField setup global view test", "[Datatypes::ScalarFieldSetup_global_view_uniformUp]")
{
   std::size_t dim3D;
   std::size_t dim1D;
   std::vector<std::size_t> idx3D;
   auto spSetup = details::createSetup(dim1D, dim3D, idx3D, details::SetupType::UniformUp, 1);
   auto&& meta = *spSetup->viewMeta();

   auto dfct = [](const std::size_t k){return k + 2;};
   auto sze = [&](const std::size_t& a, const auto& b){return a + dfct(b);};
   auto acc = [&](const std::size_t n){return std::accumulate(idx3D.begin(), idx3D.begin()+n, std::size_t(0), sze);};

   // Count layers
   std::size_t nLayers = 0;
   for(std::size_t i = 0; i < meta.ptr2D.size() - 1; i++)
   {
      if(meta.ptr2D.at(i+1) > meta.ptr2D.at(i))
      {
         nLayers++;
      }
   }

   CHECK( nLayers == idx3D.size() );
   CHECK( *std::max_element(meta.dim1D.begin(), meta.dim1D.end()) == dim1D );
   CHECK( *meta.ptr2D.rbegin() == acc(idx3D.size()) );
   CHECK( meta.ptr2D.at(idx3D.at(1)) == acc(1) );
   CHECK( meta.dim1D.at(idx3D.at(2)) == dim1D );
   CHECK( meta.ptr2D.at(idx3D.at(3)+1) - meta.ptr2D.at(idx3D.at(3)) == dfct(idx3D.at(3)) );
   CHECK( (1 + meta.ptr2D.at(idx3D.at(0))) == 1 );
   CHECK( (1 + meta.ptr2D.at(idx3D.at(0))) == 1 );
   CHECK( (2 + meta.ptr2D.at(idx3D.at(1))) == acc(1) + 2 );
   CHECK( (2 + meta.ptr2D.at(idx3D.at(3))) == acc(3) + 2 );
}

TEST_CASE("Uniform decreasing ScalarField setup test", "[Datatypes::ScalarFieldSetup_uniformDown]")
{
   std::size_t dim3D;
   std::size_t dim1D;
   std::vector<std::size_t> idx3D;
   auto spSetup = details::createSetup(dim1D, dim3D, idx3D, details::SetupType::UniformDown);

   auto dfct = [&](const std::size_t k){return dim1D - k;};
   auto sze = [&](const std::size_t& a, const auto& b){return a + dfct(b);};
   auto acc = [&](const std::size_t n){return static_cast<int>(std::accumulate(idx3D.begin(), idx3D.begin()+n, 0, sze));};

   int d3D = static_cast<int>(idx3D.size());
   int d1D = static_cast<int>(dim1D);

   CHECK( spSetup->nBlock() == d3D );
   CHECK( spSetup->dataRows() == d1D );
   CHECK( spSetup->dataCols() == acc(idx3D.size()) );
   CHECK( spSetup->blockIdx(1) == acc(1) );
   CHECK( spSetup->blockRows(2) == d1D );
   CHECK( spSetup->blockCols(3) == static_cast<int>(dfct(3)) );
   CHECK( spSetup->colIdx(3) == 3 );
   CHECK( spSetup->colIdx(3, 0) == 3 );
   CHECK( spSetup->colIdx(2, 3) == acc(3) + 2 );
}

TEST_CASE("Uniform decreasing ScalarField setup view test", "[Datatypes::ScalarFieldSetup_view_uniformDown]")
{
   std::size_t dim3D;
   std::size_t dim1D;
   std::vector<std::size_t> idx3D;
   auto spSetup = details::createSetup(dim1D, dim3D, idx3D, details::SetupType::UniformDown);
   auto&& meta = *spSetup->viewMeta();

   auto dfct = [&](const std::size_t k){return dim1D - k;};
   auto sze = [&](const std::size_t& a, const auto& b){return a + dfct(b);};
   auto acc = [&](const std::size_t n){return std::accumulate(idx3D.begin(), idx3D.begin()+n, std::size_t(0), sze);};

   CHECK( meta.ptr2D.size() - 1 == dim3D );
   CHECK( *std::max_element(meta.dim1D.begin(), meta.dim1D.end()) == dim1D );
   CHECK( *meta.ptr2D.rbegin() == acc(idx3D.size()) );
   CHECK( meta.ptr2D.at(1) == acc(1) );
   CHECK( meta.dim1D.at(2) == dim1D );
   CHECK( meta.ptr2D.at(4) - meta.ptr2D.at(3) == dfct(idx3D.at(3)) );
   CHECK( (3 + meta.ptr2D.at(0)) == 3 );
   CHECK( (3 + meta.ptr2D.at(0)) == 3 );
   CHECK( (2 + meta.ptr2D.at(3)) == acc(3) + 2 );
}

TEST_CASE("Triangular increasing ScalarField setup test", "[Datatypes::ScalarFieldSetup_triangularUp]")
{
   std::size_t dim3D;
   std::size_t dim1D;
   std::vector<std::size_t> idx3D;
   auto spSetup = details::createSetup(dim1D, dim3D, idx3D, details::SetupType::TriangularUp);

   auto dfct = [&](const std::size_t k){return k + 1;};
   auto sze = [&](const std::size_t& a, const auto& b){return a + dfct(b);};
   auto acc = [&](const std::size_t n){return static_cast<int>(std::accumulate(idx3D.begin(), idx3D.begin()+n, 0, sze));};

   int d3D = static_cast<int>(idx3D.size());
   int d1D = static_cast<int>(dim1D);

   CHECK( spSetup->nBlock() == d3D );
   CHECK( spSetup->dataRows() == d1D + d3D - 1 );
   CHECK( spSetup->dataCols() == acc(idx3D.size()) );
   CHECK( spSetup->blockIdx(1) == acc(1) );
   CHECK( spSetup->blockRows(2) == (2 + d1D) );
   CHECK( spSetup->blockCols(3) == dfct(3) );
   CHECK( spSetup->colIdx(0) == 0 );
   CHECK( spSetup->colIdx(0, 0) == 0 );
   CHECK( spSetup->colIdx(1, 1) == acc(1) + 1 );
   CHECK( spSetup->colIdx(2, 3) == acc(3) + 2 );
}

TEST_CASE("Triangular increasing ScalarField setup view test", "[Datatypes::ScalarFieldSetup_view_triangularUp]")
{
   std::size_t dim3D;
   std::size_t dim1D;
   std::vector<std::size_t> idx3D;
   auto spSetup = details::createSetup(dim1D, dim3D, idx3D, details::SetupType::TriangularUp);
   auto&& meta = *spSetup->viewMeta();

   auto dfct = [&](const std::size_t k){return k + 1;};
   auto sze = [&](const std::size_t& a, const auto& b){return a + dfct(b);};
   auto acc = [&](const std::size_t n){return std::accumulate(idx3D.begin(), idx3D.begin()+n, std::size_t(0), sze);};

   CHECK( meta.ptr2D.size() - 1 == dim3D );
   CHECK( *std::max_element(meta.dim1D.begin(), meta.dim1D.end()) == dim1D + dim3D - 1 );
   CHECK( *meta.ptr2D.rbegin() == acc(idx3D.size()) );
   CHECK( meta.ptr2D.at(1) == acc(1) );
   CHECK( meta.dim1D.at(meta.ptr2D.at(2)) == 2 + dim1D );
   CHECK( meta.ptr2D.at(4) - meta.ptr2D.at(3) == dfct(3) );
   CHECK( (0 + meta.ptr2D.at(0)) == 0 );
   CHECK( (0 + meta.ptr2D.at(0)) == 0 );
   CHECK( (1 + meta.ptr2D.at(1)) == acc(1) + 1 );
   CHECK( (2 + meta.ptr2D.at(3)) == acc(3) + 2 );
}

TEST_CASE("Triangular decreasing ScalarField setup test", "[Datatypes::ScalarFieldSetup_triangularDown]")
{
   std::size_t dim3D;
   std::size_t dim1D;
   std::vector<std::size_t> idx3D;
   auto spSetup = details::createSetup(dim1D, dim3D, idx3D, details::SetupType::TriangularDown);

   auto dfct = [&](const std::size_t k){return dim1D - 2*k;};
   auto sze = [&](const std::size_t& a, const auto& b){return a + dfct(b);};
   auto acc = [&](const std::size_t n){return static_cast<int>(std::accumulate(idx3D.begin(), idx3D.begin()+n, 0, sze));};

   int d3D = static_cast<int>(idx3D.size());
   int d1D = static_cast<int>(dim1D);

   CHECK( spSetup->nBlock() == d3D );
   CHECK( spSetup->dataRows() == d1D );
   CHECK( spSetup->dataCols() == acc(idx3D.size()) );
   CHECK( spSetup->blockIdx(1) == acc(1) );
   CHECK( spSetup->blockRows(2) == d1D - 2 );
   CHECK( spSetup->blockCols(3) == static_cast<int>(dfct(3)) );
   CHECK( spSetup->colIdx(3) == 3 );
   CHECK( spSetup->colIdx(3, 0) == 3 );
   CHECK( spSetup->colIdx(2, 3) == acc(3) + 2 );
}

TEST_CASE("Triangular decreasing ScalarField setup view test", "[Datatypes::ScalarFieldSetup_view_triangularDown]")
{
   std::size_t dim3D;
   std::size_t dim1D;
   std::vector<std::size_t> idx3D;
   auto spSetup = details::createSetup(dim1D, dim3D, idx3D, details::SetupType::TriangularDown);
   auto&& meta = *spSetup->viewMeta();

   auto dfct = [&](const std::size_t k){return dim1D - 2*k;};
   auto sze = [&](const std::size_t& a, const auto& b){return a + dfct(b);};
   auto acc = [&](const std::size_t n){return std::accumulate(idx3D.begin(), idx3D.begin()+n, std::size_t(0), sze);};

   CHECK( meta.ptr2D.size() - 1 == dim3D );
   CHECK( *std::max_element(meta.dim1D.begin(), meta.dim1D.end()) == dim1D );
   CHECK( *meta.ptr2D.rbegin() == acc(idx3D.size()) );
   CHECK( meta.ptr2D.at(1) == acc(1) );
   CHECK( meta.dim1D.at(meta.ptr2D.at(2)) == dim1D - 2 );
   CHECK( meta.ptr2D.at(4) - meta.ptr2D.at(3) == dfct(3) );
   CHECK( (3 + meta.ptr2D.at(0)) == 3 );
   CHECK( (3 + meta.ptr2D.at(0)) == 3 );
   CHECK( (2 + meta.ptr2D.at(3)) == acc(3) + 2 );
}

