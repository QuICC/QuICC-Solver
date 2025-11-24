#include <catch2/catch.hpp>

#include "QuICC/Enums/FieldIds.hpp"
#include "details/Helper.hpp"

TEST_CASE("FlatScalarField get point data", "[FlatScalarField::point]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;
   
   CHECK( field.point(0,0,0) == details::fieldValueA<double>(0,0,0) );
   CHECK( field.point(3,0,0) == details::fieldValueA<double>(3,0,0) );
   CHECK( field.point(1,1,0) == details::fieldValueA<double>(1,1,0) );
   CHECK( field.point(1,0,1) == details::fieldValueA<double>(1,0,1) );
   CHECK( field.point(3,2,3) == details::fieldValueA<double>(3,2,3) );
}

TEST_CASE("FlatScalarField get point data from coord", "[FlatScalarField::point_coord]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;
   
   std::vector<int> coord = {0,0,0};
   CHECK( field.point(coord) == details::fieldValueA<double>(0,0,0) );
   coord = {3, 0, 0};
   CHECK( field.point(coord) == details::fieldValueA<double>(3,0,0) );
   coord = {1, 1, 0};
   CHECK( field.point(coord) == details::fieldValueA<double>(1,1,0) );
   coord = {1, 0, 1};
   CHECK( field.point(coord) == details::fieldValueA<double>(1,0,1) );
   coord = {3, 2, 3};
   CHECK( field.point(coord) == details::fieldValueA<double>(3,2,3) );
}

TEST_CASE("FlatScalarField get profile data", "[FlatScalarField::profile]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;
  
   auto&& p0 = field.profile(1);
   for(int i = 0; i < dim1D; i++)
   {
      CHECK( p0(i) == details::fieldValueA<double>(i,1,0) );
   }
  
   auto&& p3 = field.profile(2, 3);
   for(int i = 0; i < dim1D; i++)
   {
      CHECK( p3(i) == details::fieldValueA<double>(i,2,3) );
   }
  
   auto&& p4 = field.profile(5, 4);
   for(int i = 0; i < dim1D; i++)
   {
      CHECK( p4(i) == details::fieldValueA<double>(i,5,4) );
   }
}

TEST_CASE("FlatScalarField get slice data", "[FlatScalarField::slice]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;
  
   int k = 0;
   auto&& s0 = field.slice(k);
   for(int j = 0; j < k + 2; j++)
   {
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( s0(i, j) == details::fieldValueA<double>(i,j,k) );
      }
   }
  
   k = 2;
   auto&& s2 = field.slice(k);
   for(int j = 0; j < k + 2; j++)
   {
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( s2(i, j) == details::fieldValueA<double>(i,j,k) );
      }
   }
  
   k = 4;
   auto&& s4 = field.slice(k);
   for(int j = 0; j < k + 2; j++)
   {
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( s4(i, j) == details::fieldValueA<double>(i,j,k) );
      }
   }
}

TEST_CASE("FlatScalarField set point data", "[FlatScalarField::setPoint]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   field.setPoint(-700.007, 0, 1, 0);
   field.setPoint(-700.008, 2, 1, 1);
   field.setPoint(-700.009, 7, 3, 3);
   field.setPoint(-700.017, dim1D-1, 1, 4);

   CHECK( field.point(0, 1, 0) == -700.007);
   CHECK( field.point(2, 1, 1) == -700.008);
   CHECK( field.point(7, 3, 3) == -700.009);
   CHECK( field.point(dim1D-1, 1, 4) == -700.017);
}

TEST_CASE("FlatScalarField set point data from variant", "[FlatScalarField::setPoint_variant]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;
   QuICC::MHDVariant pt;

   pt = -700.007;
   field.setPoint(pt, 0, 1, 0);
   pt = -700.008;
   field.setPoint(pt, 2, 1, 1);
   pt = -700.009;
   field.setPoint(pt, 7, 3, 3);
   pt = -700.017;
   field.setPoint(pt, dim1D-1, 1, 4);

   CHECK( field.point(0, 1, 0) == -700.007);
   CHECK( field.point(2, 1, 1) == -700.008);
   CHECK( field.point(7, 3, 3) == -700.009);
   CHECK( field.point(dim1D-1, 1, 4) == -700.017);
}

TEST_CASE("FlatScalarField set profile data", "[FlatScalarField::setProfile]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   QuICC::Array np(dim1D);
   for(int i = 0; i < dim1D; i++)
   {
      np(i) = details::fieldValueB<double>(i,0,0);
   }

   field.setProfile(np, 0, 0);
   field.setProfile(np, 3, 2);
   field.setProfile(np, 3, 3);
   field.setProfile(np, 1, 4);
  
   {
      int j = 0;
      int k = 0;
      auto&& p = field.profile(j, k);
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( p(i) == details::fieldValueB<double>(i,0,0) );
      }
   }
  
   {
      int j = 3;
      int k = 2;
      auto&& p = field.profile(j, k);
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( p(i) == details::fieldValueB<double>(i,0,0) );
      }
   }
  
   {
      int j = 3;
      int k = 3;
      auto&& p = field.profile(j, k);
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( p(i) == details::fieldValueB<double>(i,0,0) );
      }
   }
  
   {
      int j = 1;
      int k = 4;
      auto&& p = field.profile(j, k);
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( p(i) == details::fieldValueB<double>(i,0,0) );
      }
   }
  
   {
      int j = 2;
      int k = 2;
      auto&& p = field.profile(j, k);
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( p(i) == details::fieldValueA<double>(i,j,k) );
      }
   }
}

TEST_CASE("FlatScalarField add profile data", "[FlatScalarField::addProfile]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   QuICC::Array np(dim1D);
   for(int i = 0; i < dim1D; i++)
   {
      np(i) = details::fieldValueB<double>(i,0,0);
   }

   field.addProfile(np, 0, 0);
   field.addProfile(np, 3, 2);
   field.addProfile(np, 3, 3);
   field.addProfile(np, 1, 4);
  
   {
      int j = 0;
      int k = 0;
      auto&& p = field.profile(j, k);
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( p(i) == details::fieldValueA<double>(i,j,k) +  (details::fieldValueB<double>(i,0,0)) );
      }
   }
  
   {
      int j = 3;
      int k = 2;
      auto&& p = field.profile(j, k);
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( p(i) == details::fieldValueA<double>(i,j,k) + (details::fieldValueB<double>(i,0,0)) );
      }
   }
  
   {
      int j = 3;
      int k = 3;
      auto&& p = field.profile(j, k);
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( p(i) == details::fieldValueA<double>(i,j,k) + (details::fieldValueB<double>(i,0,0)) );
      }
   }
  
   {
      int j = 1;
      int k = 4;
      auto&& p = field.profile(j, k);
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( p(i) == details::fieldValueA<double>(i,j,k) + (details::fieldValueB<double>(i,0,0)) );
      }
   }
  
   {
      int j = 2;
      int k = 2;
      auto&& p = field.profile(j, k);
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( p(i) == details::fieldValueA<double>(i,j,k) );
      }
   }
}

TEST_CASE("FlatScalarField subtract profile data", "[FlatScalarField::subProfile]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   QuICC::Array np(dim1D);
   for(int i = 0; i < dim1D; i++)
   {
      np(i) = details::fieldValueB<double>(i,0,0);
   }

   field.subProfile(np, 0, 0);
   field.subProfile(np, 3, 2);
   field.subProfile(np, 3, 3);
   field.subProfile(np, 1, 4);
  
   {
      int j = 0;
      int k = 0;
      auto&& p = field.profile(j, k);
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( p(i) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,0,0)) );
      }
   }
  
   {
      int j = 3;
      int k = 2;
      auto&& p = field.profile(j, k);
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( p(i) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,0,0)) );
      }
   }
  
   {
      int j = 3;
      int k = 3;
      auto&& p = field.profile(j, k);
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( p(i) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,0,0)) );
      }
   }
  
   {
      int j = 1;
      int k = 4;
      auto&& p = field.profile(j, k);
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( p(i) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,0,0)) );
      }
   }
  
   {
      int j = 2;
      int k = 2;
      auto&& p = field.profile(j, k);
      for(int i = 0; i < dim1D; i++)
      {
         CHECK( p(i) == details::fieldValueA<double>(i,j,k) );
      }
   }
}

TEST_CASE("FlatScalarField set slice data", "[FlatScalarField::setSlice]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   {
      int k = 0;
      QuICC::Matrix ns(dim1D, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.setSlice(ns, k);

      k = 2;
      ns.resize(dim1D, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.setSlice(ns, k);

      k = 4;
      ns.resize(dim1D, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.setSlice(ns, k);
   }

   {
      int k = 0;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
         }
      }
   }

   {
      int k = 1;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }

   {
      int k = 2;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
         }
      }
   }

   {
      int k = 4;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
         }
      }
   }
}

TEST_CASE("FlatScalarField add slice data", "[FlatScalarField::addSlice]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   {
      int k = 0;
      QuICC::Matrix ns(dim1D, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.addSlice(ns, k);

      k = 2;
      ns.resize(dim1D, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.addSlice(ns, k);

      k = 4;
      ns.resize(dim1D, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.addSlice(ns, k);
   }

   {
      int k = 0;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) + (details::fieldValueB<double>(i,j,k)) );
         }
      }
   }

   {
      int k = 1;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }

   {
      int k = 2;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) + (details::fieldValueB<double>(i,j,k)) );
         }
      }
   }

   {
      int k = 4;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) + (details::fieldValueB<double>(i,j,k)) );
         }
      }
   }
}

TEST_CASE("FlatScalarField subtract slice data", "[FlatScalarField::subSlice]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   {
      int k = 0;
      QuICC::Matrix ns(dim1D, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.subSlice(ns, k);

      k = 2;
      ns.resize(dim1D, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.subSlice(ns, k);

      k = 4;
      ns.resize(dim1D, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.subSlice(ns, k);
   }

   {
      int k = 0;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,j,k)) );
         }
      }
   }

   {
      int k = 1;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }

   {
      int k = 2;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,j,k)) );
         }
      }
   }

   {
      int k = 4;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,j,k)) );
         }
      }
   }
}

TEST_CASE("FlatScalarField set top rows of slice data", "[FlatScalarField::setSlice_rows]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   {
      int k = 0;
      QuICC::Matrix ns(dim1D/2, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.setSlice(ns, k, dim1D/2);

      k = 2;
      ns.resize(dim1D/2, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.setSlice(ns, k, dim1D/2);

      k = 4;
      ns.resize(dim1D/2, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.setSlice(ns, k, dim1D/2);
   }

   {
      int k = 0;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
         }
         for(int i = dim1D/2; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }

   {
      int k = 1;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }

   {
      int k = 2;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
         }
         for(int i = dim1D/2; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }

   {
      int k = 4;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
         }
         for(int i = dim1D/2; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }
}

TEST_CASE("FlatScalarField add to top rows of slice data", "[FlatScalarField::addSlice_rows]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   {
      int k = 0;
      QuICC::Matrix ns(dim1D/2, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.addSlice(ns, k, dim1D/2);

      k = 2;
      ns.resize(dim1D/2, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.addSlice(ns, k, dim1D/2);

      k = 4;
      ns.resize(dim1D/2, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.addSlice(ns, k, dim1D/2);
   }

   {
      int k = 0;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) + (details::fieldValueB<double>(i,j,k)) );
         }
         for(int i = dim1D/2; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }

   {
      int k = 1;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }

   {
      int k = 2;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) + (details::fieldValueB<double>(i,j,k)) );
         }
         for(int i = dim1D/2; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }

   {
      int k = 4;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) + (details::fieldValueB<double>(i,j,k)) );
         }
         for(int i = dim1D/2; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }
}

TEST_CASE("FlatScalarField subtract from top rows of slice data", "[FlatScalarField::subSlice_rows]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   {
      int k = 0;
      QuICC::Matrix ns(dim1D/2, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.subSlice(ns, k, dim1D/2);

      k = 2;
      ns.resize(dim1D/2, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.subSlice(ns, k, dim1D/2);

      k = 4;
      ns.resize(dim1D/2, k + 2);
      ns.setZero();
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            ns(i, j) = details::fieldValueB<double>(i,j,k);
         }
      }

      field.subSlice(ns, k, dim1D/2);
   }

   {
      int k = 0;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,j,k)) );
         }
         for(int i = dim1D/2; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }

   {
      int k = 1;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }

   {
      int k = 2;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,j,k)) );
         }
         for(int i = dim1D/2; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }

   {
      int k = 4;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D/2; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,j,k)) );
         }
         for(int i = dim1D/2; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }
}

TEST_CASE("FlatScalarField get data pointer", "[FlatScalarField::data]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   {
      int k = 0;
      auto&& dptr = field.data(k);
      
      CHECK( dptr[0] == 0 );
      CHECK( dptr[3] == 3 );
      CHECK( dptr[1*dim1D + 3] == 1*100 + 3 );
   }

   {
      int k = 2;
      auto&& dptr = field.data(k);
      
      CHECK( dptr[0] == k*100*100 + 0 );
      CHECK( dptr[3] == k*100*100 + 3 );
      CHECK( dptr[1*dim1D + 1] == k*100*100 + 101 );
   }

   {
      int k = 4;
      auto&& dptr = field.data(k);
      
      CHECK( dptr[0] == k*100*100 + 0 );
      CHECK( dptr[3] == k*100*100 + 3 );
      CHECK( dptr[1*dim1D + 2] == k*100*100 + 102 );
   }
}

TEST_CASE("FlatScalarField get data reference", "[FlatScalarField::data_ref]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   {
      auto&& mat = field.data();
      
      CHECK( mat(0, 0) == 0 );
      CHECK( mat(3, 0) == 3 );
      CHECK( mat(2, 1) == 1*100 + 2 );
   }

   {
      int k = 2;
      auto&& mat = field.data();
      
      CHECK( mat(0, (0 + 2) + (1 + 2) + 0) == k*100*100 + 0 );
      CHECK( mat(3, (0 + 2) + (1 + 2) + 0) == k*100*100 + 3 );
      CHECK( mat(1, (0 + 2) + (1 + 2) + 1) == k*100*100 + 1*100 + 1 );
   }

   {
      int k = 4;
      auto&& mat = field.data();
      
      CHECK( mat(0, (0 + 2) + (1 + 2) + (2 + 2) + (3 + 2) + 0) == k*100*100 + 0 );
      CHECK( mat(3, (0 + 2) + (1 + 2) + (2 + 2) + (3 + 2) + 0) == k*100*100 + 3 );
      CHECK( mat(3, (0 + 2) + (1 + 2) + (2 + 2) + (3 + 2) + 1) == k*100*100 + 1*100 + 3 );
   }
}

TEST_CASE("FlatScalarField set data", "[FlatScalarField::setData]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   {
      QuICC::Matrix nd(dim1D, ((dim3D)*(dim3D + 3))/2);
      nd.setZero();
      int j_ = 0;
      for(int k = 0; k < dim3D; k++)
      {
         for(int j = 0; j < k + 2; j++)
         {
            for(int i = 0; i < dim1D; i++)
            {
               nd(i, j_) = details::fieldValueB<double>(i,j,k);
            }
            j_++;
         }
      }

      field.setData(nd);
   }

   {
      int k = 0;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
         }
      }
   }

   {
      int k = 2;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
         }
      }
   }

   {
      int k = 4;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
         }
      }
   }
}

TEST_CASE("FlatScalarField set data with flipped sign", "[FlatScalarField::setNegData]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   {
      QuICC::Matrix nd(dim1D, ((dim3D)*(dim3D + 3))/2);
      nd.setZero();
      int j_ = 0;
      for(int k = 0; k < dim3D; k++)
      {
         for(int j = 0; j < k + 2; j++)
         {
            for(int i = 0; i < dim1D; i++)
            {
               nd(i, j_) = details::fieldValueB<double>(i,j,k);
            }
            j_++;
         }
      }

      field.setNegData(nd);
   }

   {
      int k = 0;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == -(details::fieldValueB<double>(i,j,k)) );
         }
      }
   }

   {
      int k = 2;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == -(details::fieldValueB<double>(i,j,k)) );
         }
      }
   }

   {
      int k = 4;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == -(details::fieldValueB<double>(i,j,k)) );
         }
      }
   }
}

TEST_CASE("FlatScalarField add data", "[FlatScalarField::addData]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   {
      QuICC::Matrix nd(dim1D, ((dim3D)*(dim3D + 3))/2);
      nd.setZero();
      int j_ = 0;
      for(int k = 0; k < dim3D; k++)
      {
         for(int j = 0; j < k + 2; j++)
         {
            for(int i = 0; i < dim1D; i++)
            {
               nd(i, j_) = details::fieldValueB<double>(i,j,k);
            }
            j_++;
         }
      }

      field.addData(nd);
   }

   {
      int k = 0;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) + (details::fieldValueB<double>(i,j,k)) );
         }
      }
   }

   {
      int k = 2;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) + (details::fieldValueB<double>(i,j,k)) );
         }
      }
   }

   {
      int k = 4;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) + (details::fieldValueB<double>(i,j,k)) );
         }
      }
   }
}

TEST_CASE("FlatScalarField subtract data", "[FlatScalarField::subData]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   {
      QuICC::Matrix nd(dim1D, ((dim3D)*(dim3D + 3))/2);
      nd.setZero();
      int j_ = 0;
      for(int k = 0; k < dim3D; k++)
      {
         for(int j = 0; j < k + 2; j++)
         {
            for(int i = 0; i < dim1D; i++)
            {
               nd(i, j_) = details::fieldValueB<double>(i,j,k);
            }
            j_++;
         }
      }

      field.subData(nd);
   }

   {
      int k = 0;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,j,k)) );
         }
      }
   }

   {
      int k = 2;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,j,k)) );
         }
      }
   }

   {
      int k = 4;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,j,k)) );
         }
      }
   }
}

TEST_CASE("FlatScalarField set zeros", "[FlatScalarField::setZeros]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   field.setZeros();

   {
      int k = 0;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == 0 );
         }
      }
   }

   {
      int k = 2;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == 0 );
         }
      }
   }

   {
      int k = 4;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == 0 );
         }
      }
   }
}

TEST_CASE("FlatScalarField rescale data by scalar", "[FlatScalarField::rescale]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   double c = -3.1415;
   field.rescale(c);

   std::vector<int> ks = {0, 2, 4};
   for(auto&& k: ks)
   {
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == c*(details::fieldValueA<double>(i,j,k)) );
         }
      }
   }
}

TEST_CASE("FlatScalarField get number of slices", "[FlatScalarField::nSlice]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   CHECK( field.nSlice() == dim3D );
}

TEST_CASE("FlatScalarField set slice data from pointer", "[FlatScalarField::rData]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   {
      int k = 0;
      double* dptr = field.rData(k);

      int ii = 0;
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            dptr[ii] = details::fieldValueB<double>(i,j,k);
            ii++;
         }
      }

      k = 2;
      dptr = field.rData(k);
      ii = 0;
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            dptr[ii] = details::fieldValueB<double>(i,j,k);
            ii++;
         }
      }

      k = 4;
      dptr = field.rData(k);
      ii = 0;
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            dptr[ii] = details::fieldValueB<double>(i,j,k);
            ii++;
         }
      }
   }

   {
      int k = 0;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
         }
      }
   }

   {
      int k = 1;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }

   {
      int k = 2;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
         }
      }
   }

   {
      int k = 4;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
         }
      }
   }
}

TEST_CASE("FlatScalarField set whole data from reference", "[FlatScalarField::rData_ref]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   {
      QuICC::Matrix nd(dim1D, ((dim3D + 1)*(dim3D + 4))/2);
      nd.setZero();
      int j_ = 0;
      for(int k = 0; k < dim3D; k++)
      {
         for(int j = 0; j < k + 2; j++)
         {
            for(int i = 0; i < dim1D; i++)
            {
               nd(i, j_) = details::fieldValueB<double>(i,j,k);
            }
            j_++;
         }
      }

      field.rData() = nd;
   }

   {
      int k = 0;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
         }
      }
   }

   {
      int k = 2;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
         }
      }
   }

   {
      int k = 4;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(int i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
         }
      }
   }
}

TEST_CASE("FlatScalarField set point data from reference", "[FlatScalarField::rPoint]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   field.rPoint(0, 1, 0) = -700.007;
   field.rPoint(2, 1, 1) = -700.008;
   field.rPoint(7, 3, 3) = -700.009;
   field.rPoint(dim1D-1, 1, 4) = -700.017;

   CHECK( field.point(0, 1, 0) == -700.007);
   CHECK( field.point(2, 1, 1) == -700.008);
   CHECK( field.point(7, 3, 3) == -700.009);
   CHECK( field.point(dim1D-1, 1, 4) == -700.017);
}

TEST_CASE("FlatScalarField set point data from reference from coord", "[FlatScalarField::rPoint_coord]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   std::vector<int> coord = {0, 1, 0};
   field.rPoint(coord) = -700.007;
   coord = {2, 1, 1};
   field.rPoint(coord) = -700.008;
   coord = {7, 3, 3};
   field.rPoint(coord) = -700.009;
   coord = {dim1D-1, 1, 4};
   field.rPoint(coord) = -700.017;

   CHECK( field.point(0, 1, 0) == -700.007);
   CHECK( field.point(2, 1, 1) == -700.008);
   CHECK( field.point(7, 3, 3) == -700.009);
   CHECK( field.point(dim1D-1, 1, 4) == -700.017);
}

TEST_CASE("FlatScalarField copy constructor", "[FlatScalarField::Copy]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spFieldA = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);
   auto spFieldB = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);
   auto spFieldC = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& fieldA = *spFieldA;
   fieldA.setZeros();
   {
      auto ptr = fieldA.data().data();

      for(int i = 0; i < fieldA.data().size(); i++)
      {
         CHECK( *ptr == 0 );
         ptr++;
      }
   }

   auto&& fieldB = *spFieldB;
   auto&& fieldC = *spFieldC;

   fieldA = fieldB;
   {
      auto ptrA = fieldA.data().data();
      auto ptrC = fieldC.data().data();

      for(int i = 0; i < fieldA.data().size(); i++)
      {
         CHECK( *ptrA == *ptrC );
         ptrA++;
         ptrC++;
      }
   }

   fieldB.setZeros();
   {
      auto ptrA = fieldA.data().data();

      for(int i = 0; i < fieldA.data().size(); i++)
      {
         CHECK( *ptrA == 0 );
         ptrA++;
      }
   }
}

TEST_CASE("FlatScalarField get component by reference", "[FlatScalarField::comp]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   CHECK( &(field.comp(QuICC::FieldComponents::Spectral::SCALAR)) == &field );
}

TEST_CASE("FlatScalarField set component by reference", "[FlatScalarField::rComp]")
{
   int dim3D = 5;
   int dim1D = 2*dim3D;
   auto spField = details::createFlatScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D);

   auto&& field = *spField;

   CHECK( &(field.rComp(QuICC::FieldComponents::Spectral::SCALAR)) == &field );
}
