/**
 * @file SphericalCoriolis.hpp
 * @brief Implementation of the spherical coriolis term
 */

#ifndef QUICC_PHYSICAL_SPHERICALCORIOLIS_HPP
#define QUICC_PHYSICAL_SPHERICALCORIOLIS_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "ViewOps/Slicewise/NoGridOp.hpp"
#include "ViewOps/Slicewise/NoTwoGridOp.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of the spherical coriolis term
    */
   class SphericalCoriolis
   {
      public:
         /**
          * @brief Set S to Coriolis term
          */
         template <typename TFIELD>
            static void set(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis term to S
          */
         template <typename TFIELD>
            static void add(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Substract Coriolis term from S
          */
         template <typename TFIELD>
            static void sub(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Set S to Coriolis term
          */
         template <typename TFIELD, typename TIDXFUNC>
            static void set(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis term to S
          */
         template <typename TFIELD, typename TIDXFUNC>
            static void add(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Substract Coriolis term from S
          */
         template <typename TFIELD, typename TIDXFUNC>
            static void sub(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalCoriolis() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalCoriolis() = default;

      private:
         /// Functor to map resolution object
         struct IdxResFunctor
         {
            const Resolution& _res;

            const int dim3D;

            IdxResFunctor(const Resolution& res) : _res(res), dim3D(res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>()) {};

            /// @brief deleted default constructor
            IdxResFunctor() = delete;

            /// @brief dtor
            ~IdxResFunctor() = default;

            int dim2D(const int k) const
            {
               return _res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(k);
            }

            int idx2D(const int j, const int k) const
            {
               return _res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(j, k);
            }
         };

         /// @tparam T scalar
         template <class T = double> struct SetRTFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            SetRTFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            SetRTFunctor() = delete;

            /// @brief dtor
            ~SetRTFunctor() = default;

            /// @brief Dot product
            /// @param g
            /// @param ui
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T g, T ui)
            {
               return -_scaling * (g * ui);
            }
         };

         /// @tparam T scalar
         template <class T = double> struct AddRTFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            AddRTFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            AddRTFunctor() = delete;

            /// @brief dtor
            ~AddRTFunctor() = default;

            /// @brief Dot product
            /// @param g
            /// @param ui
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T g, T w, T ui)
            {
               return w - _scaling * (g * ui);
            }
         };

         /// @tparam T scalar
         template <class T = double> struct SubRTFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            SubRTFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            SubRTFunctor() = delete;

            /// @brief dtor
            ~SubRTFunctor() = default;

            /// @brief Dot product
            /// @param g
            /// @param ui
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T g, T w, T ui)
            {
               return w + _scaling * (g * ui);
            }
         };

         /// @tparam T scalar
         template <class T = double> struct SetPFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            SetPFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            SetPFunctor() = delete;

            /// @brief dtor
            ~SetPFunctor() = default;

            /// @brief Dot product
            /// @param gc
            /// @param gs
            /// @param ui
            /// @param uj
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T gc, T gs, T ui, T uj)
            {
               return _scaling * (gc * ui + gs * uj);
            }
         };

         /// @tparam T scalar
         template <class T = double> struct AddPFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            AddPFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            AddPFunctor() = delete;

            /// @brief dtor
            ~AddPFunctor() = default;

            /// @brief Dot product
            /// @param gc
            /// @param gs
            /// @param ui
            /// @param uj
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T gc, T gs, T w, T ui, T uj)
            {
               return w + _scaling * (gc * ui + gs * uj);
            }
         };

         /// @tparam T scalar
         template <class T = double> struct SubPFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            SubPFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            SubPFunctor() = delete;

            /// @brief dtor
            ~SubPFunctor() = default;

            /// @brief Dot product
            /// @param gc
            /// @param gs
            /// @param ui
            /// @param uj
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T gc, T gs, T w, T ui, T uj)
            {
               return w - _scaling * (gc * ui + gs * uj);
            }
         };
   };

   template <typename TFIELD>
   void SphericalCoriolis::set(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      IdxResFunctor f(res);
      set(rS, compId, f, cosTheta, sinTheta, v, c);
   }

   template <typename TFIELD>
   void SphericalCoriolis::add(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      IdxResFunctor f(res);
      add(rS, compId, f, cosTheta, sinTheta, v, c);
   }

   template <typename TFIELD>
   void SphericalCoriolis::sub(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      IdxResFunctor f(res);
      sub(rS, compId, f, cosTheta, sinTheta, v, c);
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalCoriolis::set(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      int nR = idxFunc.dim3D;
      int nTh;
      int iTh_;

      if(compId == FieldComponents::Physical::R)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = SetRTFunctor<scalar_t>;
            fct_t f(c);
            grid_t vSin(const_cast<scalar_t *>(sinTheta.data()), sinTheta.size());
            Slicewise::Cpu::NoGridOp<1, fct_t, view_t, grid_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vSin, v.comp(FieldComponents::Physical::PHI).globalView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.setProfile(-c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.setProfile(-v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
                  }
               }
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = SetRTFunctor<scalar_t>;
            fct_t f(c);
            grid_t vCos(const_cast<scalar_t *>(cosTheta.data()), cosTheta.size());
            Slicewise::Cpu::NoGridOp<1, fct_t, view_t, grid_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vCos, v.comp(FieldComponents::Physical::PHI).globalView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.setProfile(-c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.setProfile(-v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
                  }
               }
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = SetPFunctor<scalar_t>;
            fct_t f(c);
            grid_t vCos(const_cast<scalar_t *>(cosTheta.data()), cosTheta.size());
            grid_t vSin(const_cast<scalar_t *>(sinTheta.data()), sinTheta.size());
            Slicewise::Cpu::NoTwoGridOp<1, fct_t, view_t, grid_t, grid_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vCos, vSin, v.comp(FieldComponents::Physical::THETA).globalView(), v.comp(FieldComponents::Physical::R).globalView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.setProfile(c*(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
                     rS.addProfile(c*(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.setProfile(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
                     rS.addProfile(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
                  }
               }
            }
         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalCoriolis::add(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      int nR = idxFunc.dim3D;
      int nTh;
      int iTh_;

      if(compId == FieldComponents::Physical::R)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = AddRTFunctor<scalar_t>;
            fct_t f(c);
            grid_t vSin(const_cast<scalar_t *>(sinTheta.data()), sinTheta.size());
            Slicewise::Cpu::NoGridOp<1, fct_t, view_t, grid_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vSin, rS.dataView(), v.comp(FieldComponents::Physical::PHI).globalView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.subProfile(c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.subProfile(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
                  }
               }
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = AddRTFunctor<scalar_t>;
            fct_t f(c);
            grid_t vCos(const_cast<scalar_t *>(cosTheta.data()), cosTheta.size());
            Slicewise::Cpu::NoGridOp<1, fct_t, view_t, grid_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vCos, rS.dataView(), v.comp(FieldComponents::Physical::PHI).globalView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.subProfile(c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.subProfile(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
                  }
               }
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = AddPFunctor<scalar_t>;
            fct_t f(c);
            grid_t vCos(const_cast<scalar_t *>(cosTheta.data()), cosTheta.size());
            grid_t vSin(const_cast<scalar_t *>(sinTheta.data()), sinTheta.size());
            Slicewise::Cpu::NoTwoGridOp<1, fct_t, view_t, grid_t, grid_t, view_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vCos, vSin, rS.dataView(), v.comp(FieldComponents::Physical::THETA).globalView(), v.comp(FieldComponents::Physical::R).globalView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.addProfile(c*(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
                     rS.addProfile(c*(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.addProfile(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
                     rS.addProfile(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
                  }
               }
            }
         }
      }
   }
   
   template <typename TFIELD, typename TIDXFUNC>
      void SphericalCoriolis::sub(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      int nR = idxFunc.dim3D;
      int nTh;
      int iTh_;

      if(compId == FieldComponents::Physical::R)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = SubRTFunctor<scalar_t>;
            fct_t f(c);
            grid_t vSin(const_cast<scalar_t *>(sinTheta.data()), sinTheta.size());
            Slicewise::Cpu::NoGridOp<1, fct_t, view_t, grid_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vSin, rS.dataView(), v.comp(FieldComponents::Physical::PHI).globalView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.addProfile(c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.addProfile(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
                  }
               }
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = SubRTFunctor<scalar_t>;
            fct_t f(c);
            grid_t vCos(const_cast<scalar_t *>(cosTheta.data()), cosTheta.size());
            Slicewise::Cpu::NoGridOp<1, fct_t, view_t, grid_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vCos, rS.dataView(), v.comp(FieldComponents::Physical::PHI).globalView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.addProfile(c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.addProfile(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
                  }
               }
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = SubPFunctor<scalar_t>;
            fct_t f(c);
            grid_t vCos(const_cast<scalar_t *>(cosTheta.data()), cosTheta.size());
            grid_t vSin(const_cast<scalar_t *>(sinTheta.data()), sinTheta.size());
            Slicewise::Cpu::NoTwoGridOp<1, fct_t, view_t, grid_t, grid_t, view_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vCos, vSin, rS.dataView(), v.comp(FieldComponents::Physical::THETA).globalView(), v.comp(FieldComponents::Physical::R).globalView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.subProfile(c*(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
                     rS.subProfile(c*(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  nTh = idxFunc.dim2D(iR); 
                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.subProfile(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
                     rS.subProfile(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
                  }
               }
            }
         }
      }
   }
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_SPHERICALCORIOLIS_HPP
