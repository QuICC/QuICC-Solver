#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

#include "details/Helper.hpp"

TEST_CASE("Uniform increasing ScalarField setup test", "[Datatypes::ScalarFieldSetup_uniformUp]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spSetup = details::createSetup(details::SetupType::UniformUp, dim1D, dim3D);

   CHECK( spSetup->nBlock() == dim3D );
   CHECK( spSetup->dataRows() == dim1D );
   CHECK( spSetup->dataCols() == (0 + 2) + (1 + 2) + (2 + 2) + (3 + 2) + (4 + 2) );
   CHECK( spSetup->blockIdx(1) == (0 + 2) );
   CHECK( spSetup->blockRows(2) == dim1D );
   CHECK( spSetup->blockCols(3) == 3 + 2 );
   CHECK( spSetup->colIdx(3) == 3 );
   CHECK( spSetup->colIdx(3, 0) == 3 );
   CHECK( spSetup->colIdx(2, 3) == (0 + 2) + (1 + 2) + (2 + 2) + 2 );
}

TEST_CASE("Uniform decreasing ScalarField setup test", "[Datatypes::ScalarFieldSetup_uniformDown]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spSetup = details::createSetup(details::SetupType::UniformDown, dim1D, dim3D);

   CHECK( spSetup->nBlock() == dim3D );
   CHECK( spSetup->dataRows() == dim1D );
   CHECK( spSetup->dataCols() == (dim1D - 0) + (dim1D - 1) + (dim1D - 2) + (dim1D - 3) + (dim1D - 4) );
   CHECK( spSetup->blockIdx(1) == (dim1D - 0) );
   CHECK( spSetup->blockRows(2) == dim1D );
   CHECK( spSetup->blockCols(3) == dim1D - 3 );
   CHECK( spSetup->colIdx(3) == 3 );
   CHECK( spSetup->colIdx(3, 0) == 3 );
   CHECK( spSetup->colIdx(2, 3) == (dim1D - 0) + (dim1D - 1) + (dim1D - 2) + 2 );
}

TEST_CASE("Triangular increasing ScalarField setup test", "[Datatypes::ScalarFieldSetup_triangularUp]")
{
   int dim3D = 5;
   int dim1D = 2;
   auto spSetup = details::createSetup(details::SetupType::TriangularUp, dim1D, dim3D);

   CHECK( spSetup->nBlock() == dim3D );
   CHECK( spSetup->dataRows() == dim1D + dim3D - 1 );
   CHECK( spSetup->dataCols() == (0 + 1) + (1 + 1) + (2 + 1) + (3 + 1) + (4 + 1) );
   CHECK( spSetup->blockIdx(1) == (0 + 1) );
   CHECK( spSetup->blockRows(2) == (2 + dim1D) );
   CHECK( spSetup->blockCols(3) == 3 + 1 );
   CHECK( spSetup->colIdx(3) == 3 );
   CHECK( spSetup->colIdx(3, 0) == 3 );
   CHECK( spSetup->colIdx(2, 3) == (0 + 1) + (1 + 1) + (2 + 1) + 2 );
}

TEST_CASE("Triangular decreasing ScalarField setup test", "[Datatypes::ScalarFieldSetup_triangularDown]")
{
   int dim3D = 5;
   int dim1D = 3*dim3D;
   auto spSetup = details::createSetup(details::SetupType::TriangularDown, dim1D, dim3D);

   CHECK( spSetup->nBlock() == dim3D );
   CHECK( spSetup->dataRows() == dim1D );
   CHECK( spSetup->dataCols() == (dim1D - 2*0) + (dim1D - 2*1) + (dim1D - 2*2) + (dim1D - 2*3) + (dim1D - 2*4) );
   CHECK( spSetup->blockIdx(1) == (dim1D - 2*0) );
   CHECK( spSetup->blockRows(2) == dim1D - 2 );
   CHECK( spSetup->blockCols(3) == dim1D - 2*3 );
   CHECK( spSetup->colIdx(3) == 3 );
   CHECK( spSetup->colIdx(3, 0) == 3 );
   CHECK( spSetup->colIdx(2, 3) == (dim1D - 2*0) + (dim1D - 2*1) + (dim1D - 2*2) + 2 );
}

