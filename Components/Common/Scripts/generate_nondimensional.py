#!/usr/bin/env python

def header(descr):
    cName = descr['class']
    name = descr['name']
    content = f'''/**
 * @file {cName}.hpp
 * @brief {name} number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_{cName.upper()}_HPP
#define QUICC_NONDIMENSIONAL_{cName.upper()}_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/BasicTypes.hpp"
#include "QuICC/NonDimensional/IRegisterId.hpp"

namespace QuICC {{

namespace NonDimensional {{

   /**
    * @brief {name} number nondimensional number
    */
   class {cName}: public IRegisterId<{cName}>
   {{
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of {name} number
          */
         {cName}(const MHDFloat value);

         friend class IRegisterId<{cName}>;

      protected:

      private:
         /**
          * @brief Unique tag
          */
         static std::string sTag();

         /**
          * @brief Formatted name
          */
         static std::string sFormatted();
   }};

}}
}}

#endif // QUICC_NONDIMENSIONAL_{cName.upper()}_HPP
'''
    return content

def source(descr):
    cName = descr['class']
    name = descr['name']
    tag = descr['tag']
    content = f'''/**
 * @file {cName}.cpp
 * @brief Source of the {name} nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/{cName}.hpp"

// Project includes
//

namespace QuICC {{

namespace NonDimensional {{

   std::string {cName}::sTag()
   {{
      return "{tag}";
   }}

   std::string {cName}::sFormatted()
   {{
      return "{name}";
   }}

   {cName}::{cName}(const MHDFloat value)
      : IRegisterId<{cName}>(value, {cName}::sTag(), {cName}::sFormatted())
   {{
   }}

}}
}}
'''
    return content

names = []
names.append({'class':"Alpha", 'name':"Alpha", 'tag':"alpha"})
names.append({'class':"Beta", 'name':"Beta", 'tag':"beta"})
names.append({'class':"Alpha", 'name':"Alpha", 'tag':"alpha"})
names.append({'class':"CflAlfvenDamping", 'name':"Alfven damping for CFL", 'tag':"cfl_alfven_damping"})
names.append({'class':"CflAlfvenScale", 'name':"Alfven scale for CFL", 'tag':"cfl_alfven_scale"})
names.append({'class':"CflInertial", 'name':"Inertial CFL", 'tag':"cfl_inertial"})
names.append({'class':"CflTorsional", 'name':"Torsional CFL", 'tag':"cfl_torsional"})
names.append({'class':"Chandrasekhar", 'name':"Chandrasekhar", 'tag':"chandrasekhar"})
names.append({'class':"Chi", 'name':"Chi", 'tag':"chi"})
names.append({'class':"Delta", 'name':"Delta", 'tag':"delta"})
names.append({'class':"Eady", 'name':"Eady", 'tag':"eady"})
names.append({'class':"Ekman", 'name':"Ekman", 'tag':"ekman"})
names.append({'class':"Elevator", 'name':"Elevator", 'tag':"elevator"})
names.append({'class':"Elsasser", 'name':"Elsasser", 'tag':"elsasser"})
names.append({'class':"Epsilon", 'name':"Epsilon", 'tag':"epsilon"})
names.append({'class':"Eta", 'name':"Eta", 'tag':"eta"})
names.append({'class':"FastMean", 'name':"fast mean", 'tag':"fast_mean"})
names.append({'class':"Gamma", 'name':"Gamma", 'tag':"gamma"})
names.append({'class':"Heating", 'name':"Heating", 'tag':"heating"})
names.append({'class':"Iota", 'name':"Iota", 'tag':"iota"})
names.append({'class':"Kappa", 'name':"Kappa", 'tag':"kappa"})
names.append({'class':"Lambda", 'name':"Lambda", 'tag':"lambda"})
names.append({'class':"Lower1D", 'name':"Lower1D", 'tag':"lower1d"})
names.append({'class':"Lower2D", 'name':"Lower2D", 'tag':"lower2d"})
names.append({'class':"Lower3D", 'name':"Lower3D", 'tag':"lower3d"})
names.append({'class':"MagEkman", 'name':"magnetic Ekman", 'tag':"magnetic_ekman"})
names.append({'class':"MagPrandtl", 'name':"magnetic Prandtl", 'tag':"magnetic_prandtl"})
names.append({'class':"MagReynolds", 'name':"magnetic Reynolds", 'tag':"magnetic_reynolds"})
names.append({'class':"ModElsasser", 'name':"modified Elsasser", 'tag':"modified_elsasser"})
names.append({'class':"Mu", 'name':"Mu", 'tag':"mu"})
names.append({'class':"Nu", 'name':"Nu", 'tag':"nu"})
names.append({'class':"Omega", 'name':"Omega", 'tag':"omega"})
names.append({'class':"Omicron", 'name':"Omicron", 'tag':"omicron"})
names.append({'class':"Phi", 'name':"Phi", 'tag':"phi"})
names.append({'class':"Pi", 'name':"Pi", 'tag':"pi"})
names.append({'class':"Poincare", 'name':"Poincare", 'tag':"poincare"})
names.append({'class':"Prandtl", 'name':"Prandtl", 'tag':"prandtl"})
names.append({'class':"Psi", 'name':"Psi", 'tag':"psi"})
names.append({'class':"RRatio", 'name':"R ratio", 'tag':"rratio"})
names.append({'class':"Rayleigh", 'name':"Rayleigh", 'tag':"rayleigh"})
names.append({'class':"Rescaled", 'name':"Rescaled", 'tag':"rescaled"})
names.append({'class':"Rho", 'name':"Rho", 'tag':"rho"})
names.append({'class':"Roberts", 'name':"Roberts", 'tag':"roberts"})
names.append({'class':"Rossby", 'name':"Rossby", 'tag':"rossby"})
names.append({'class':"Sigma", 'name':"Sigma", 'tag':"sigma"})
names.append({'class':"Tau", 'name':"Tau", 'tag':"tau"})
names.append({'class':"Taylor", 'name':"Taylor", 'tag':"taylor"})
names.append({'class':"Theta", 'name':"Theta", 'tag':"theta"})
names.append({'class':"Upper1D", 'name':"Upper1D", 'tag':"upper1d"})
names.append({'class':"Upper2D", 'name':"Upper2D", 'tag':"upper2d"})
names.append({'class':"Upper3D", 'name':"Upper3D", 'tag':"upper3d"})
names.append({'class':"Upsilon", 'name':"Upsilon", 'tag':"upsilon"})
names.append({'class':"Xi", 'name':"Xi", 'tag':"xi"})
names.append({'class':"Zeta", 'name':"Zeta", 'tag':"zeta"})

incPath = '../include/QuICC/NonDimensional/'
srcPath = '../src/NonDimensional/'
for n in names:
    f = open(incPath+n['class']+'.hpp', 'w')
    f.write(header(n))
    f.close()
    f = open(srcPath+n['class']+'.cpp', 'w')
    f.write(source(n))
    f.close()
