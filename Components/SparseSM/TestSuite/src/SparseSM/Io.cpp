/** 
 * @file Io.cpp
 * @brief Source of the useful functions for SparseSM tests
 */

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <unsupported/Eigen/SparseExtra>

// External includes
//

// Project includes
//
#include "QuICC/Math/Constants.hpp"
#include "QuICC/TestSuite/SparseSM/Io.hpp"

namespace QuICC {

namespace TestSuite {

namespace SparseSM {

   void writeData(const std::string& path, const SparseMatrix& outData)
   {
      Eigen::saveMarket(outData, path);
   }

   void writeData(const std::string& path, const Matrix& outData)
   {
      std::ofstream outfile;
      outfile.open(path);
      if(! outfile.is_open())
      {
         throw std::logic_error("Couldn't open output file!");
      }
      outfile << std::setprecision(15) << outData << std::endl;
      outfile.close();
   }

   void readData(SparseMatrix& inData, const std::string& path)
   {
      // Check if file exists
      std::ifstream infile;
      infile.open(path, std::ios::in | std::ios::binary);
      if(! infile.is_open())
      {
         throw std::logic_error("Couldn't open input file: " + path);
      }
      infile.close();

      Eigen::loadMarket(inData, path);
   }

   void readData(Matrix& inData, const std::string& path)
   {
      std::ifstream infile;
      infile.open(path, std::ios::in | std::ios::binary);
      if(! infile.is_open())
      {
         std::cerr << "*****************************************************************" << std::endl;
         std::cerr << "*****************************************************************" << std::endl;
         std::cerr << "  Couldn't open real input file: " + path << std::endl;
         std::cerr << "*****************************************************************" << std::endl;
         std::cerr << "*****************************************************************" << std::endl;
         inData.resize(0,0);
      } else
      {
         // Get size from first two values
         if(inData.size() == 0)
         {
            MHDFloat r;
            infile >> r;
            MHDFloat c;
            infile >> c;
            inData.resize(static_cast<int>(r), static_cast<int>(c));
         }

         // Loop over data
         for(int i = 0; i < inData.rows(); ++i)
         {
            for(int j = 0; j < inData.cols(); ++j)
            {
               infile >> inData(i,j);
            }
         }
         infile.close();
      }
   }

   void readList(Array& inData, const std::string& path)
   {
      std::ifstream infile;
      infile.open(path, std::ios::in | std::ios::binary);
      if(! infile.is_open())
      {
         std::cerr << "*****************************************************************" << std::endl;
         std::cerr << "*****************************************************************" << std::endl;
         std::cerr << "  Couldn't open real input file: " + path << std::endl;
         std::cerr << "*****************************************************************" << std::endl;
         std::cerr << "*****************************************************************" << std::endl;
         inData.setConstant(std::numeric_limits<MHDFloat>::max());
      } else
      {
         // Get size from first value
         if(inData.size() == 0)
         {
            MHDFloat s;
            infile >> s;
            inData.resize(static_cast<int>(s));
         }

         // Loop over data
         for(int i = 0; i < inData.size(); ++i)
         {
            infile >> inData(i);
         }
         infile.close();
      }
   }

}
}
}
