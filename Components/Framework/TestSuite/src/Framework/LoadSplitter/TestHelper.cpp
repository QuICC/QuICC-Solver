/**
 * @file TestHelper.cpp
 * @brief Source of test helper
 */

// System includes
//
#include <catch2/catch.hpp>
#include <iostream>
#include <fstream>
#include <limits>
#include <set>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/TestSuite/Framework/LoadSplitter/TestHelper.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace LoadSplitter {

   void writeData(const std::string fname, const SimRes& sRes, const QuICC::TransformResolution& tRes)
   {
      int nModes = tRes.dim<QuICC::Dimensions::Data::DAT3D>();

      std::size_t nIndexes = 0;
      for(auto k = 0; k < tRes.dim<QuICC::Dimensions::Data::DAT3D>(); k++)
      {
         nIndexes += tRes.dim<QuICC::Dimensions::Data::DAT2D>(k);
      }

      std::ofstream metaFile;
      metaFile.open(fname);
      if(! metaFile.is_open())
      {
         throw std::logic_error("Couldn't open output file!");
      }

      // Write header
      metaFile << "# Format description" << std::endl;
      metaFile << "# line 1: number of lines" << std::endl;
      metaFile << "# line 2: simulation modal size" << std::endl;
      metaFile << "# line 3: simulation physical size" << std::endl;
      metaFile << "# line 4: number of modes nK" << std::endl;
      metaFile << "# line 5 + 2*k: mode" << std::endl;
      metaFile << "# line 6 + 2*k: number of 2D modes" << std::endl;
      metaFile << "# line (5 + 2*nK) + 2*j: 2D mode" << std::endl;
      metaFile << "# line (6 + 2*nK) + 2*j: number of 1D modes" << std::endl;

      metaFile << 3 + 2*nModes + 2*nIndexes << std::endl;
      if(nModes > 0)
      {
         metaFile << sRes.nSpec << std::endl;
         metaFile << sRes.nPhys << std::endl;
         metaFile << nModes << std::endl;
      }
      else
      {
         metaFile << 0 << std::endl;;
         metaFile << 0 << std::endl;;
         metaFile << 0 << std::endl;;
      }

      for(int i = 0; i < nModes; i++)
      {
         int m = tRes.idx<QuICC::Dimensions::Data::DAT3D>(i);
         int mult = tRes.dim<QuICC::Dimensions::Data::DAT2D>(i);

         if(metaFile.is_open())
         {
            metaFile << m << std::endl;
            metaFile << mult << std::endl;
         }
      }

      for(auto k = 0; k < tRes.dim<QuICC::Dimensions::Data::DAT3D>(); k++)
      {
         for(int j = 0; j < tRes.dim<QuICC::Dimensions::Data::DAT2D>(k); j++)
         {
            auto j_ = tRes.idx<QuICC::Dimensions::Data::DAT2D>(j, k);
            auto nN = tRes.dim<QuICC::Dimensions::Data::DATB1D>(j, k);
            metaFile << j_ << std::endl;
            metaFile << nN << std::endl;
         }
      }

      if(metaFile.is_open())
      {
         metaFile.close();
      }
   }

   void readData(const std::string path, std::vector<int>& data)
   {
      data.clear();

      std::ifstream infile;
      infile.open(path, std::ios::in | std::ios::binary);
      if(! infile.is_open())
      {
         std::cerr << "*****************************************************************" << std::endl;
         std::cerr << "*****************************************************************" << std::endl;
         std::cerr << "  Couldn't open input metadata file: " + path << std::endl;
         std::cerr << "*****************************************************************" << std::endl;
         std::cerr << "*****************************************************************" << std::endl;
      }
      else
      {
         // Ignore header
         int s = infile.peek();
         while(s == '#')
         {
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            s = infile.peek();
         }

         // Read size and prepare storage
         infile >> s;
         data.reserve(s);

         // Loop over data
         int v;
         for(int i = 0; i < s; ++i)
         {
            infile >> v;
            data.push_back(v);
         }
         infile.close();
      }
   }

   void collectModes(const QuICC::TransformResolution& tRes, std::multimap<int,int>& modes)
   {
      // Collect (i,j) modes
      for(int k = 0; k < tRes.dim<QuICC::Dimensions::Data::DAT3D>(); k++)
      {
         auto k_ = tRes.idx<QuICC::Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < tRes.dim<QuICC::Dimensions::Data::DAT2D>(k); j++)
         {
            auto j_ = tRes.idx<QuICC::Dimensions::Data::DAT2D>(j, k);
            modes.emplace(k_, j_);
         }
      }
   }

   void checkReference(const std::string fname, const SimRes& sRes, const QuICC::TransformResolution& tRes)
   {
      std::vector<int> ref;
      readData(fname, ref);

      // Check reference data exists
      CHECK( ref.size() > 0 );

      if(ref.size() > 0)
      {
         // Check 1D sizes
         {
            INFO( "Checking 1D arrays sizes" );
            if(tRes.dim<QuICC::Dimensions::Data::DAT3D>() > 0)
            {
               // Definition of size to use is changing (ie including dealiasing modes or not)
               CHECK( ref.at(0) == sRes.nSpec );
               CHECK( ref.at(1) == sRes.nPhys );
            }
            else
            {
               CHECK( ref.at(0) == 0 );
               CHECK( ref.at(1) == 0 );
            }
         }

         // Check number of 3D modes
         {
            INFO( "Checking number of 3D modes" );
            CHECK( ref.at(2) == tRes.dim<QuICC::Dimensions::Data::DAT3D>() );
         }

         // Check number of 2D modes per 3D mode
         {
            INFO( "Checking number of 2D modes per 3D index" );
            for(int k = 0; k < tRes.dim<QuICC::Dimensions::Data::DAT3D>(); k++)
            {
               auto jump = k*2;
               auto k_ = tRes.idx<QuICC::Dimensions::Data::DAT3D>(k);
               auto nJ = tRes.dim<QuICC::Dimensions::Data::DAT2D>(k);
               auto ref_k_ = ref.at(3 + jump);
               auto ref_nJ = ref.at(4 + jump);
               INFO( "ref k = " << ref_k_ );
               INFO( "k = " << k_ );
               CHECK( ref_k_ == k_ );
               INFO( "ref nJ = " << ref_nJ );
               INFO( "nJ = " << nJ );
               CHECK( ref_nJ == nJ );
            }
         }

         // Check details of modes
         auto h = 2 + 2*tRes.dim<QuICC::Dimensions::Data::DAT3D>() + 1;
         for(int k = 0; k < tRes.dim<QuICC::Dimensions::Data::DAT3D>(); k++)
         {
            for(int j = 0; j < tRes.dim<QuICC::Dimensions::Data::DAT2D>(k); j++)
            {
               auto j_ = tRes.idx<QuICC::Dimensions::Data::DAT2D>(j, k);
               auto nI = tRes.dim<QuICC::Dimensions::Data::DATB1D>(j, k);
               if(static_cast<std::size_t>(h) < ref.size())
               {
                  auto ref_j_ = ref.at(h);
                  auto ref_nI = ref.at(h+1);
                  INFO( "ref j = " << ref_j_ );
                  INFO( "j = " << j_ );
                  CHECK( ref_j_ == j_ );
                  INFO( "ref nI = " << ref_nI );
                  INFO( "nI = " << nI );
                  CHECK( ref_nI == nI );
               }
               else
               {
                  INFO( "Point missing in reference" );
                  CHECK( false );
               }
               h += 2;
            }
         }
      }
   }

   void checkSerialReference(const std::string fname, const std::multimap<int,int>& modes)
   {
      std::vector<int> ref;
      readData(fname, ref);

      // Check reference data exists
      CHECK( ref.size() > 0 );

      if(ref.size() > 0)
      {
         // Number of 3D modes
         auto n3D = ref.at(2);

         // Check number of 2D modes per 3D mode
         {
            INFO( "Checking number of 2D modes per 3D index against Serial" );
            for(int k = 0; k < n3D; k++)
            {
               auto jump = k*2;
               auto ref_k_ = ref.at(3 + jump);
               std::size_t ref_nJ = ref.at(4 + jump);
               auto nJ = modes.count(ref_k_);
               INFO( "ref k = " << ref_k_ );
               INFO( "ref nJ = " << ref_nJ );
               INFO( "nJ = " << nJ );
               CHECK( ref_nJ == nJ );
            }
         }

         // Check details of modes
         {
            INFO( "Checking 2D modes per 3D index against Serial" );
            auto h = 2 + 2*n3D + 1;
            for(int k = 0; k < n3D; k++)
            {
               auto jump = k*2;
               auto k_ = ref.at(3 + jump);
               auto r = modes.equal_range(k_);
               std::set<int> js;
               for(auto it = r.first; it != r.second; ++it)
               {
                  js.insert(it->second);
               }
               // Checking duplicate modes
               {
                  INFO( "check for duplicate modes" );
                  CHECK( js.size() == modes.count(k_) );
               }

               for(auto j_: js)
               {
                  auto ref_j_ = ref.at(h);
                  INFO( "ref j = " << ref_j_ );
                  INFO( "j = " << j_ );
                  CHECK( ref_j_ == j_ );
                  h += 2;
               }
            }
         }
      }
   }

} // LoadSplitter
} // Framework
} // TestSuite
} // QuICC
