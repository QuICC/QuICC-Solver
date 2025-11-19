#include <catch2/catch.hpp>

#include "QuICC/Variables/Spectral/ScalarVariable.hpp"
#include "details/Helper.hpp"

TEST_CASE("FlatScalarField get point data", "[FlatScalarField::point]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;
   
   CHECK( field.point(0,0,0) == 0 );
   CHECK( field.point(3,0,0) == 3 );
   CHECK( field.point(1,1,0) == 101 );
   CHECK( field.point(1,0,1) == 10001 );
   CHECK( field.point(3,2,3) == 30203 );
}

TEST_CASE("FlatScalarField get point data from coord", "[FlatScalarField::point_coord]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;
   
   std::vector<int> coord = {0,0,0};
   CHECK( field.point(coord) == 0 );
   coord = {3, 0, 0};
   CHECK( field.point(coord) == 3 );
   coord = {1, 1, 0};
   CHECK( field.point(coord) == 101 );
   coord = {1, 0, 1};
   CHECK( field.point(coord) == 10001 );
   coord = {3, 2, 3};
   CHECK( field.point(coord) == 30203 );
}

TEST_CASE("FlatScalarField get profile data", "[FlatScalarField::profile]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField get slice data", "[FlatScalarField::slice]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField set point data", "[FlatScalarField::setPoint]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField set point data from variant", "[FlatScalarField::setPoint_variant]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField set profile data", "[FlatScalarField::setProfile]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField add profile data", "[FlatScalarField::addProfile]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField subtract profile data", "[FlatScalarField::subProfile]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField set slice data", "[FlatScalarField::setSlice]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField multiply slice data by scalar", "[FlatScalarField::multSlice]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField add slice data", "[FlatScalarField::addSlice]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField subtract slice data", "[FlatScalarField::subSlice]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField set top rows of slice data", "[FlatScalarField::setSlice_rows]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField multiply top rows of slice data by scalar", "[FlatScalarField::multSlice_rows]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField add to top rows of slice data", "[FlatScalarField::addSlice_rows]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField subtract from top rows of slice data", "[FlatScalarField::subSlice_rows]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField get data pointer", "[FlatScalarField::data]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField get data reference", "[FlatScalarField::data_ref]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField set data", "[FlatScalarField::setData]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField add data", "[FlatScalarField::addData]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField subtract data", "[FlatScalarField::subData]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField set zeros", "[FlatScalarField::setZeros]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField rescale data by scalar", "[FlatScalarField::rescale]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField get number of slices", "[FlatScalarField::nSlice]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField set whole data from pointer", "[FlatScalarField::rData]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField set whole data from reference", "[FlatScalarField::rData_ref]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField set point data from reference", "[FlatScalarField::rPoint]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField set point data from reference of coord", "[FlatScalarField::rPoint_coord]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField get component by reference", "[FlatScalarField::comp]")
{
   const int a = 2;
}

TEST_CASE("FlatScalarField set component by reference", "[FlatScalarField::rComp]")
{
   const int a = 2;
}
