/**
 * @file HighamExponential.cpp
 * @brief Implementation of exponential matrix algorithm from higham
 */

// System includes
//
#include <Eigen/Dense>

// Project includes
//
#include "Timestep/Exponential/HighamExponential.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

extern "C" void dgebal_(char* job, int* n, double* A, int* lda, int* ilo,
   int* ihi, double* scale, int* info);

extern "C" void dgebak_(char* job, char* side, int* n, int* ilo, int* ihi,
   double* scale, int* m, double* A, int* lda, int* info);

extern "C" double dlange_(char* norm, int* m, int* n, double* A, int* lda,
   double* work);

HighamExponential::HighamExponential()
{
   // Set Pade coefficients
   this->initPade();
}

void HighamExponential::initPade()
{
   // Pade of order 3
   this->mPade3 = {120., 60., 12., 1.};

   // Pade of order 5
   this->mPade5 = {30240., 15120., 3360., 420., 30., 1.};

   // Pade of order 7
   this->mPade7 = {17297280., 8648640., 1995840., 277200., 25200., 1512., 56.,
      1.};

   // Pade of order 9
   this->mPade9 = {17643225600., 8821612800., 2075673600., 302702400.,
      30270240., 2162160., 110880., 3960., 90., 1.};

   // Pade of order 13
   this->mPade13 = {64764752532480000., 32382376266240000., 7771770303897600.,
      1187353796428800., 129060195264000., 10559470521600., 670442572800.,
      33522128640., 1323241920., 40840800., 960960., 16380., 182., 1.};

   // Theta_m thresholds
   this->mTheta = {1.495585217958292e-2, 2.539398330063230e-1,
      9.504178996162932e-1, 2.097847961257068e0, 5.371920351148152e0};
}

Matrix HighamExponential::compute(const Matrix& matA) const
{
   assert(matA.rows() == matA.cols());
   int n = matA.rows();

   Matrix expA = matA;

   // preprocess A to reduce norm
   std::vector<double> scale;
   int ilo, ihi;
   this->preprocess(scale, ilo, ihi, expA);

   auto normA = this->norm(matA);

   // Use Pade-approximation
   if (normA < this->mTheta[3])
   {
      const double* b;
      int bn;
      // Pade order 3
      if (normA < this->mTheta[0])
      {
         b = this->mPade3.data();
         bn = this->mPade3.size();
      }
      // Pade order 5
      else if (normA < this->mTheta[1])
      {
         b = this->mPade5.data();
         bn = this->mPade5.size();
      }
      // Pade order 7
      else if (normA < this->mTheta[2])
      {
         b = this->mPade7.data();
         bn = this->mPade7.size();
      }
      // Pade order 9
      else
      {
         b = this->mPade9.data();
         bn = this->mPade9.size();
      }

      Matrix matA2 = expA * expA;
      Matrix matQ = Matrix::Zero(n, n);
      Matrix matU = Matrix::Zero(n, n);
      Matrix matV = Matrix::Zero(n, n);

      for (int i = 0; i < n; i++)
      {
         matU(i, i) = b[1];
         matV(i, i) = b[0];
      }

      int K = bn / 2 - 1;
      for (int k = 1, k2 = 2; k <= K; k++, k2 = 2 * k)
      {
         if (k == 1)
         {
            matQ = matA2;
         }
         else
         {
            matQ = matQ * matA2;
         }

         matU += b[k2 + 1] * matQ;
         matV += b[k2] * matQ;
      }

      matU = expA * matU;
      matQ = matV - matU;
      matU += matV;

      expA = matQ.fullPivLu().solve(matU);
   }
   // Norm is too large, use scaling and squaring
   else
   {
      auto&& b = this->mPade13;

      // Scale matrix to reduce norm
      int si = std::ceil(std::log2(normA / 5.4));
      if (si > 0)
      {
         expA.array() /= std::pow(2.0, si);
      }

      Matrix matA2 = expA * expA;
      Matrix matA4 = matA2 * matA2;
      Matrix matA6 = matA2 * matA4;

      Matrix matI = Matrix::Identity(n, n);

      Matrix matU =
         expA * (matA6 * (b[13] * matA6 + b[11] * matA4 + b[9] * matA2) +
                   b[7] * matA6 + b[5] * matA4 + b[3] * matA2 + b[1] * matI);
      Matrix matV = matA6 * (b[12] * matA6 + b[10] * matA4 + b[8] * matA2) +
                    b[6] * matA6 + b[4] * matA4 + b[2] * matA2 + b[0] * matI;

      Matrix matQ = matV - matU;
      matU += matV;

      expA = matQ.fullPivLu().solve(matU);
      if (si > 0)
      {
         for (int t = 1; t <= si; t++)
         {
            expA = expA * expA;
         }
      }
   }

   // Undo preprocessing
   this->postprocess(expA, scale, ilo, ihi);

   return expA;
}

void HighamExponential::preprocess(std::vector<double>& scale, int& ilo,
   int& ihi, Matrix& matA) const
{
   assert(matA.rows() == matA.cols());

   int n = matA.rows();
   int info;
   char job = 'B'; // B=Both scale and permute
   scale.resize(n);
   dgebal_(&job, &n, matA.data(), &n, &ilo, &ihi, &scale[0], &info);

   if (info != 0)
   {
      throw std::logic_error("Balancing of matrix failed");
   }

   // Rescale for zero indexing
   ilo--;
   ihi--;
   for(int i = 0; i < ilo; i++)
   {
      scale[i]--;
   }
   for(int i = ihi+1; i < n; i++)
   {
      scale[i]--;
   }
}

void HighamExponential::postprocess(Matrix& matA, std::vector<double>& scale,
   int ilo, int ihi) const
{
   assert(matA.rows() == matA.cols());

   int n = matA.rows();

   // undo scaling
   double s;
   for(int j = ilo; j <= ihi; j++)
   {
      s = scale.at(j);
      for(int i = 0; i < n; i++)
      {
         if(i == j)
         {
            continue;
         }
         matA(j,i) *= s;
         matA(i,j) /= s;
      }
   }

   // Undo the permutations
   int jj;
   if(ilo > 0)
   {
      for(int j = ilo-1; j >= 0; j--)
      {
         jj = static_cast<int>(scale.at(j));
         for(int k = 0; k < n; k++)
         {
            std::swap(matA(k,j), matA(k,jj));
         }
         for(int k = 0; k < n; k++)
         {
            std::swap(matA(j,k), matA(jj,k));
         }
      }
   }

   if(ihi < n-1)
   {
      for(int j = ihi+1; j < n; j++)
      {
         jj = static_cast<int>(scale.at(j));
         for(int k = 0; k < n; k++)
         {
            std::swap(matA(k,j), matA(k,jj));
         }
         for(int k = 0; k < n; k++)
         {
            std::swap(matA(j,k), matA(jj,k));
         }
      }
   }
}

double HighamExponential::norm(const Matrix& matA) const
{
   assert(matA.rows() == matA.cols());

   int n = matA.rows();
   char job = '1'; // Compute 1-norm

   double work;
   double norm =
      dlange_(&job, &n, &n, const_cast<double*>(matA.data()), &n, &work);
   return norm;
}

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
