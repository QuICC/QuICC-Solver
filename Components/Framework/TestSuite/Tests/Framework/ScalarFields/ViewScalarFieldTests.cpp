#include <catch2/catch.hpp>

#include "QuICC/Enums/FieldIds.hpp"
#include "details/Helper.hpp"

TEST_CASE("ViewScalarField get point data", "[ViewScalarField::point]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   std::vector<std::vector<int>> pts = {{0,0,0}, {3,0,0}, {1,1,0}, {1,0,1}, {3,2,3}};
   for(auto&& pt: pts)
   {
      CHECK( field.point(pt.at(0),pt.at(1),pt.at(2)) == details::fieldValueA<double>(pt.at(0),pt.at(1),pt.at(2)) );
   }
}

TEST_CASE("ViewScalarField get point data from coord", "[ViewScalarField::point_coord]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   std::vector<std::vector<int>> pts = {{0,0,0}, {3,0,0}, {1,1,0}, {1,0,1}, {3,2,3}};
   for(auto&& pt: pts)
   {
      CHECK( field.point(pt) == details::fieldValueA<double>(pt.at(0),pt.at(1),pt.at(2)) );
   }
}

TEST_CASE("ViewScalarField get profile data", "[ViewScalarField::profile]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   std::vector<std::pair<int,int>> jks = {{1,0},{2,3},{5,4}};
   for(auto&& jk: jks)
   {
      auto&& p = field.profile(jk.first, jk.second);
      for(std::size_t i = 0; i < dim1D; i++)
      {
         CHECK( p(i) == details::fieldValueA<double>(i,jk.first,jk.second) );
      }
   }
}

TEST_CASE("ViewScalarField get slice data", "[ViewScalarField::slice]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   std::vector<int> ks = {0, 2, 4};
   for(auto&& k: ks)
   {
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(std::size_t i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }
}

TEST_CASE("ViewScalarField set point data", "[ViewScalarField::setPoint]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

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

TEST_CASE("ViewScalarField set point data from variant", "[ViewScalarField::setPoint_variant]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

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

TEST_CASE("ViewScalarField set profile data (View)", "[ViewScalarField::setProfileView]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   std::vector<double> pData(dim1D);
   QuICC::View::ViewBase<double> pView(pData.data(), pData.size());
   std::array<std::uint32_t, 1> dims = {static_cast<std::uint32_t>(dim1D)};
   QuICC::View::View<double, QuICC::View::dense1D> np(pView, dims);
   for(std::size_t i = 0; i < dim1D; i++)
   {
      np(i) = details::fieldValueB<double>(i,0,0);
   }

   std::vector<std::pair<int,int>> jks = {{0,0}, {3,2}, {3,3}, {1,4}};
   for(auto jk: jks)
   {
      field.setProfile(np, jk.first, jk.second);
   }

   {
      std::vector<std::pair<int,int>> jks = {{0,0}, {3,2}, {3,3}, {1,4}};
      for(auto jk: jks)
      {
         auto j = jk.first;
         auto k = jk.second;
         auto&& p = field.profile(j, k);
         for(std::size_t i = 0; i < dim1D; i++)
         {
            CHECK( p(i) == details::fieldValueB<double>(i,0,0) );
         }
      }
   }

   {
      std::vector<std::pair<int,int>> jks = {{2,2}};
      for(auto jk: jks)
      {
         auto j = jk.first;
         auto k = jk.second;
         auto&& p = field.profile(j, k);
         for(std::size_t i = 0; i < dim1D; i++)
         {
            CHECK( p(i) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }
}

TEST_CASE("ViewScalarField set slice data (View)", "[ViewScalarField::setSliceView]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      std::vector<double> sdata(dim1D*(dim3D + 2));
      QuICC::View::ViewBase<double> sView(sdata.data(), sdata.size());

      std::vector<int> ks = {0, 2, 4};

      for(int k: ks)
      {
         std::array<std::uint32_t, 2> dims = {static_cast<std::uint32_t>(dim1D), static_cast<std::uint32_t>(k + 2)};
         QuICC::View::View<double, QuICC::View::dense2D> ns(sView, {dims});
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               ns(i, j) = details::fieldValueB<double>(i,j,k);
            }
         }

         field.setSlice(ns, k);
      }
   }

   {
      std::vector<int> ks = {0, 2, 4};
      for(int k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
            }
         }
      }
   }

   {
      std::vector<int> ks = {1};
      for(int k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
            }
         }
      }
   }

}

TEST_CASE("ViewScalarField set top rows of slice data (View)", "[ViewScalarField::setSliceView_rows]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      std::vector<double> sdata(dim1D*(dim3D + 2));
      QuICC::View::ViewBase<double> sView(sdata.data(), sdata.size());

      std::vector<int> ks = {0, 2, 4};

      for(int k: ks)
      {
         std::array<std::uint32_t, 2> dims = {static_cast<std::uint32_t>(dim1D/2), static_cast<std::uint32_t>(k + 2)};
         QuICC::View::View<double, QuICC::View::dense2D> ns(sView, {dims});
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D/2; i++)
            {
               ns(i, j) = details::fieldValueB<double>(i,j,k);
            }
         }

         field.setSlice(ns, k, dim1D/2);
      }
   }

   {
      std::vector<int> ks = {0,2,4};
      for(int k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D/2; i++)
            {
               CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
            }
            for(std::size_t i = dim1D/2; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
            }
         }
      }
   }

   {
      std::vector<int> ks = {1};
      for(int k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
            }
         }
      }
   }
}

TEST_CASE("ViewScalarField set data (View)", "[ViewScalarField::setDataView]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      std::vector<std::uint32_t> pointers = {0};
      std::vector<std::uint32_t> indices;

      QuICC::Matrix dMat(dim1D, (dim3D*(dim3D + 3))/2);
      dMat.setZero();
      int j_ = 0;
      for(std::size_t k = 0; k < dim3D; k++)
      {
         for(std::size_t j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               dMat(i, j_) = details::fieldValueB<double>(i,j,k);
            }
            j_++;
            indices.push_back(static_cast<std::uint32_t>(j));
         }
         pointers.push_back(static_cast<std::uint32_t>(j_));
      }

      std::array<std::uint32_t, 3> dims{static_cast<std::uint32_t>(dim1D), static_cast<std::uint32_t>(dim3D) - 1 + 2, static_cast<std::uint32_t>(dim3D)};
      QuICC::View::ViewBase<std::uint32_t> vPointers[3];
      QuICC::View::ViewBase<std::uint32_t> vIndices[3];
      vPointers[1] = QuICC::View::ViewBase<std::uint32_t>(pointers.data(), pointers.size());
      vIndices[1] = QuICC::View::ViewBase<std::uint32_t>(indices.data(), indices.size());

      QuICC::View::View<double, QuICC::View::DCCSC3D> nd(dMat.data(), dMat.size(), dims.data(), vPointers, vIndices);

      field.setData(nd);
   }

   {
      std::vector<int> ks = {0, 2, 4};
      for(auto&& k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
            }
         }
      }
   }
}

TEST_CASE("ViewScalarField set profile data (Eigen)", "[ViewScalarField::setProfile]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   QuICC::Array np(dim1D);
   for(std::size_t i = 0; i < dim1D; i++)
   {
      np(i) = details::fieldValueB<double>(i,0,0);
   }

   std::vector<std::pair<int,int>> jks = {{0,0}, {3,2}, {3,3}, {1,4}};
   for(auto jk: jks)
   {
      field.setProfile(np, jk.first, jk.second);
   }

   {
      std::vector<std::pair<int,int>> jks = {{0,0}, {3,2}, {3,3}, {1,4}};
      for(auto jk: jks)
      {
         int j = jk.first;
         int k = jk.second;
         auto&& p = field.profile(j, k);
         for(std::size_t i = 0; i < dim1D; i++)
         {
            CHECK( p(i) == (details::fieldValueB<double>(i,0,0)) );
         }
      }
   }

   {
      std::vector<std::pair<int,int>> jks = {{2,2}};
      for(auto jk: jks)
      {
         int j = jk.first;
         int k = jk.second;
         auto&& p = field.profile(j, k);
         for(std::size_t i = 0; i < dim1D; i++)
         {
            CHECK( p(i) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }
}

TEST_CASE("ViewScalarField add profile data (Eigen)", "[ViewScalarField::addProfile]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   QuICC::Array np(dim1D);
   for(std::size_t i = 0; i < dim1D; i++)
   {
      np(i) = details::fieldValueB<double>(i,0,0);
   }

   std::vector<std::pair<int,int>> jks = {{0,0}, {3,2}, {3,3}, {1,4}};
   for(auto jk: jks)
   {
      field.addProfile(np, jk.first, jk.second);
   }

   {
      std::vector<std::pair<int,int>> jks = {{0,0}, {3,2}, {3,3}, {1,4}};
      for(auto jk: jks)
      {
         int j = jk.first;
         int k = jk.second;
         auto&& p = field.profile(j, k);
         for(std::size_t i = 0; i < dim1D; i++)
         {
            CHECK( p(i) == details::fieldValueA<double>(i,j,k) +  (details::fieldValueB<double>(i,0,0)) );
         }
      }
   }

   {
      std::vector<std::pair<int,int>> jks = {{2,2}};
      for(auto jk: jks)
      {
         int j = jk.first;
         int k = jk.second;
         auto&& p = field.profile(j, k);
         for(std::size_t i = 0; i < dim1D; i++)
         {
            CHECK( p(i) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }
}

TEST_CASE("ViewScalarField subtract profile data (Eigen)", "[ViewScalarField::subProfile]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   QuICC::Array np(dim1D);
   for(std::size_t i = 0; i < dim1D; i++)
   {
      np(i) = details::fieldValueB<double>(i,0,0);
   }

   std::vector<std::pair<int,int>> jks = {{0,0}, {3,2}, {3,3}, {1,4}};
   for(auto jk: jks)
   {
      field.subProfile(np, jk.first, jk.second);
   }

   {
      std::vector<std::pair<int,int>> jks = {{0,0}, {3,2}, {3,3}, {1,4}};
      for(auto jk: jks)
      {
         int j = jk.first;
         int k = jk.second;
         auto&& p = field.profile(j, k);
         for(std::size_t i = 0; i < dim1D; i++)
         {
            CHECK( p(i) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,0,0)) );
         }
      }
   }

   {
      std::vector<std::pair<int,int>> jks = {{2,2}};
      for(auto jk: jks)
      {
         int j = jk.first;
         int k = jk.second;
         auto&& p = field.profile(j, k);
         for(std::size_t i = 0; i < dim1D; i++)
         {
            CHECK( p(i) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }
}

TEST_CASE("ViewScalarField set slice data (Eigen)", "[ViewScalarField::setSlice]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      std::vector<int> ks = {0, 2, 4};
      for(int k: ks)
      {
         QuICC::Matrix ns(dim1D, k + 2);
         ns.setZero();
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               ns(i, j) = details::fieldValueB<double>(i,j,k);
            }
         }

         field.setSlice(ns, k);
      }
   }

   {
      std::vector<int> ks = {0,2,4};
      for(int k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == (details::fieldValueB<double>(i,j,k)) );
            }
         }
      }
   }

   {
      std::vector<int> ks = {1};
      for(int k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
            }
         }
      }
   }
}

TEST_CASE("ViewScalarField add slice data (Eigen)", "[ViewScalarField::addSlice]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      std::vector<int> ks = {0, 2, 4};
      for(int k: ks)
      {
         QuICC::Matrix ns(dim1D, k + 2);
         ns.setZero();
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               ns(i, j) = details::fieldValueB<double>(i,j,k);
            }
         }

         field.addSlice(ns, k);
      }
   }

   {
      std::vector<int> ks = {0,2,4};
      for(int k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) + (details::fieldValueB<double>(i,j,k)) );
            }
         }
      }
   }

   {
      std::vector<int> ks = {1};
      for(int k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
            }
         }
      }
   }
}

TEST_CASE("ViewScalarField subtract slice data (Eigen)", "[ViewScalarField::subSlice]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      std::vector<int> ks = {0, 2, 4};
      for(int k: ks)
      {
         QuICC::Matrix ns(dim1D, k + 2);
         ns.setZero();
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               ns(i, j) = details::fieldValueB<double>(i,j,k);
            }
         }

         field.subSlice(ns, k);
      }
   }

   {
      std::vector<int> ks = {0,2,4};
      for(int k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,j,k)) );
            }
         }
      }
   }

   {
      std::vector<int> ks = {1};
      for(int k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
            }
         }
      }
   }
}

TEST_CASE("ViewScalarField set top rows of slice data (Eigen)", "[ViewScalarField::setSlice_rows]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      std::vector<int> ks = {0, 2, 4};

      for(int k: ks)
      {
         QuICC::Matrix ns(dim1D, k + 2);
         ns.setZero();
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D/2; i++)
            {
               ns(i, j) = details::fieldValueB<double>(i,j,k);
            }
         }

         field.setSlice(ns, k, dim1D/2);
      }
   }

   {
      std::vector<int> ks = {0,2,4};
      for(int k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D/2; i++)
            {
               CHECK( s(i, j) == (details::fieldValueB<double>(i,j,k)) );
            }
            for(std::size_t i = dim1D/2; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
            }
         }
      }
   }

   {
      std::vector<int> ks = {1};
      for(int k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
            }
         }
      }
   }
}

TEST_CASE("ViewScalarField add to top rows of slice data (Eigen)", "[ViewScalarField::addSlice_rows]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      std::vector<int> ks = {0, 2, 4};

      for(int k: ks)
      {
         QuICC::Matrix ns(dim1D, k + 2);
         ns.setZero();
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D/2; i++)
            {
               ns(i, j) = details::fieldValueB<double>(i,j,k);
            }
         }

         field.addSlice(ns, k, dim1D/2);
      }
   }

   {
      std::vector<int> ks = {0,2,4};
      for(int k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D/2; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) + (details::fieldValueB<double>(i,j,k)) );
            }
            for(std::size_t i = dim1D/2; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
            }
         }
      }
   }

   {
      std::vector<int> ks = {1};
      for(int k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
            }
         }
      }
   }
}

TEST_CASE("ViewScalarField subtract from top rows of slice data (Eigen)", "[ViewScalarField::subSlice_rows]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      std::vector<int> ks = {0, 2, 4};

      for(int k: ks)
      {
         QuICC::Matrix ns(dim1D, k + 2);
         ns.setZero();
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D/2; i++)
            {
               ns(i, j) = details::fieldValueB<double>(i,j,k);
            }
         }

         field.subSlice(ns, k, dim1D/2);
      }
   }

   {
      std::vector<int> ks = {0,2,4};
      for(int k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D/2; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,j,k)) );
            }
            for(std::size_t i = dim1D/2; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
            }
         }
      }
   }

   {
      int k = 1;
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(std::size_t i = 0; i < dim1D; i++)
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
         for(std::size_t i = 0; i < dim1D/2; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,j,k)) );
         }
         for(std::size_t i = dim1D/2; i < dim1D; i++)
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
         for(std::size_t i = 0; i < dim1D/2; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,j,k)) );
         }
         for(std::size_t i = dim1D/2; i < dim1D; i++)
         {
            CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
         }
      }
   }
}

TEST_CASE("ViewScalarField get data pointer", "[ViewScalarField::data]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      std::vector<std::vector<int>> ijks = {{0,0,0},{3,0,0},{2,1,0},{0,0,2},{3,0,2},{2,1,2},{0,0,4},{3,0,4},{2,1,4}};
      for(auto&& ijk: ijks)
      {
         auto&& i = ijk.at(0);
         auto&& j = ijk.at(1);
         auto&& k = ijk.at(2);
         auto&& dptr = field.data(k);

         CHECK( dptr[j*dim1D + i] == details::fieldValueA<double>(i, j, k) );
      }
   }
}

TEST_CASE("ViewScalarField get data view", "[ViewScalarField::dataView]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      auto&& v = field.dataView();

      std::vector<std::vector<int>> ijks = {{0,0,0},{3,0,0},{2,1,0},{0,0,2},{3,0,2},{2,1,2},{0,0,4},{3,0,4},{2,1,4}};
      for(auto&& ijk: ijks)
      {
         auto&& i = ijk.at(0);
         auto&& j = ijk.at(1);
         auto&& k = ijk.at(2);
         CHECK( v(i, j, k) == details::fieldValueA<double>(i, j, k) );
      }
   }
}

TEST_CASE("ViewScalarField get data (Eigen)", "[ViewScalarField::data]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      auto&& m = field.data();

      std::vector<std::vector<int>> ijks = {{0,0,0},{3,0,0},{2,1,0},{0,0,2},{3,0,2},{2,1,2},{0,0,4},{3,0,4},{2,1,4}};
      for(auto&& ijk: ijks)
      {
         auto&& i = ijk.at(0);
         auto&& j = ijk.at(1);
         auto&& k = ijk.at(2);
         auto jj = j;
         for(int kk = 0; kk < k; kk++)
         {
            jj += kk + 2;
         }
         CHECK( m(i, jj) == details::fieldValueA<double>(i, j, k) );
      }
   }
}

TEST_CASE("ViewScalarField set data (Eigen)", "[ViewScalarField::setData]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      QuICC::Matrix nd(dim1D, (dim3D*(dim3D + 3))/2);
      nd.setZero();
      int j_ = 0;
      for(std::size_t k = 0; k < dim3D; k++)
      {
         for(std::size_t j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               nd(i, j_) = details::fieldValueB<double>(i,j,k);
            }
            j_++;
         }
      }

      field.setData(nd);
   }

   {
      std::vector<int> ks = {0, 2, 4};
      for(auto&& k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
            }
         }
      }
   }
}

TEST_CASE("ViewScalarField set data with flipped sign (Eigen)", "[ViewScalarField::setNegData]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      QuICC::Matrix nd(dim1D, (dim3D*(dim3D + 3))/2);
      nd.setZero();
      int j_ = 0;
      for(std::size_t k = 0; k < dim3D; k++)
      {
         for(std::size_t j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               nd(i, j_) = details::fieldValueB<double>(i,j,k);
            }
            j_++;
         }
      }

      field.setNegData(nd);
   }

   {
      std::vector<int> ks = {0, 2, 4};
      for(auto&& k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == -(details::fieldValueB<double>(i,j,k)) );
            }
         }
      }
   }
}

TEST_CASE("ViewScalarField add data (Eigen)", "[ViewScalarField::addData]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      QuICC::Matrix nd(dim1D, (dim3D*(dim3D + 3))/2);
      nd.setZero();
      int j_ = 0;
      for(std::size_t k = 0; k < dim3D; k++)
      {
         for(std::size_t j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               nd(i, j_) = details::fieldValueB<double>(i,j,k);
            }
            j_++;
         }
      }

      field.addData(nd);
   }

   {
      std::vector<int> ks = {0, 2, 4};
      for(auto&& k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) + (details::fieldValueB<double>(i,j,k)) );
            }
         }
      }
   }
}

TEST_CASE("ViewScalarField subtract data (Eigen)", "[ViewScalarField::subData]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      QuICC::Matrix nd(dim1D, (dim3D*(dim3D + 3))/2);
      nd.setZero();
      int j_ = 0;
      for(std::size_t k = 0; k < dim3D; k++)
      {
         for(std::size_t j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               nd(i, j_) = details::fieldValueB<double>(i,j,k);
            }
            j_++;
         }
      }

      field.subData(nd);
   }

   {
      std::vector<int> ks = {0, 2, 4};
      for(auto&& k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) - (details::fieldValueB<double>(i,j,k)) );
            }
         }
      }
   }
}

TEST_CASE("ViewScalarField set zeros", "[ViewScalarField::setZeros]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   field.setZeros();

   {
      std::vector<int> ks = {0, 2, 4};
      for(auto&& k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == 0 );
            }
         }
      }
   }
}

TEST_CASE("ViewScalarField rescale data by scalar", "[ViewScalarField::rescale]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   double c = -3.1415;
   field.rescale(c);

   std::vector<int> ks = {0, 2, 4};
   for(auto&& k: ks)
   {
      auto&& s = field.slice(k);
      for(int j = 0; j < k + 2; j++)
      {
         for(std::size_t i = 0; i < dim1D; i++)
         {
            CHECK( s(i, j) == c*(details::fieldValueA<double>(i,j,k)) );
         }
      }
   }
}

TEST_CASE("ViewScalarField get number of slices", "[ViewScalarField::nSlice]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   CHECK( field.nSlice() == static_cast<int>(dim3D) );
}

TEST_CASE("ViewScalarField set slice data from pointer", "[ViewScalarField::rData]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      std::vector<int> ks = {0, 2, 4};
      for(auto&& k: ks)
      {
         double* dptr = field.rData(k);

         int ii = 0;
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               dptr[ii] = details::fieldValueB<double>(i,j,k);
               ii++;
            }
         }
      }
   }

   {
      std::vector<int> ks = {0, 2, 4};
      for(auto&& k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
            }
         }
      }
   }

   {
      std::vector<int> ks = {1};
      for(auto&& k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueA<double>(i,j,k) );
            }
         }
      }
   }
}

TEST_CASE("ViewScalarField set whole data from reference", "[ViewScalarField::rData_ref]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   {
      QuICC::Matrix nd(dim1D, ((dim3D)*(dim3D + 3))/2);
      nd.setZero();
      int j_ = 0;
      for(std::size_t k = 0; k < dim3D; k++)
      {
         for(std::size_t j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               nd(i, j_) = details::fieldValueB<double>(i,j,k);
            }
            j_++;
         }
      }

      field.rData() = nd;
   }

   {
      std::vector<int> ks = {0, 2, 4};
      for(auto&& k: ks)
      {
         auto&& s = field.slice(k);
         for(int j = 0; j < k + 2; j++)
         {
            for(std::size_t i = 0; i < dim1D; i++)
            {
               CHECK( s(i, j) == details::fieldValueB<double>(i,j,k) );
            }
         }
      }
   }
}

TEST_CASE("ViewScalarField set point data from reference", "[ViewScalarField::rPoint]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

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

TEST_CASE("ViewScalarField set point data from reference from coord", "[ViewScalarField::rPoint_coord]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   std::vector<int> coord = {0, 1, 0};
   field.rPoint(coord) = -700.007;
   coord = {2, 1, 1};
   field.rPoint(coord) = -700.008;
   coord = {7, 3, 3};
   field.rPoint(coord) = -700.009;
   coord = {static_cast<int>(dim1D-1), 1, 4};
   field.rPoint(coord) = -700.017;

   CHECK( field.point(0, 1, 0) == -700.007);
   CHECK( field.point(2, 1, 1) == -700.008);
   CHECK( field.point(7, 3, 3) == -700.009);
   CHECK( field.point(dim1D-1, 1, 4) == -700.017);
}

TEST_CASE("ViewScalarField copy constructor", "[ViewScalarField::Copy]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spFieldA = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);
   auto spFieldB = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);
   auto spFieldC = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

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

TEST_CASE("ViewScalarField get component by reference", "[ViewScalarField::comp]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   CHECK( &(field.comp(QuICC::FieldComponents::Spectral::SCALAR)) == &field );
}

TEST_CASE("ViewScalarField set component by reference", "[ViewScalarField::rComp]")
{
   std::size_t dim3D = 5;
   std::size_t dim1D = 2*dim3D;
   std::vector<std::size_t> idx3D = {0, 1, 2, 3, 4};
   auto spField = details::createViewScalarField<double>(details::SetupType::UniformUp, dim1D, dim3D, idx3D);

   auto&& field = *spField;

   CHECK( &(field.rComp(QuICC::FieldComponents::Spectral::SCALAR)) == &field );
}
