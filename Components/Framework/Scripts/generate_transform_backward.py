#!/usr/bin/env python

def header(descr):
    cName = descr['class']
    name = descr['name']
    content = f'''/**
 * @file {cName}.hpp
 * @brief Backward transform operator {name}
 */

#ifndef QUICC_TRANSFORM_BACKWARD_{cName.upper()}_HPP
#define QUICC_TRANSFORM_BACKWARD_{cName.upper()}_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Transform/Backward/IRegisterId.hpp"

namespace QuICC {{

namespace Transform {{

namespace Backward {{

   /**
    * @brief Backward transform operator {name}
    */
   class {cName}: public IRegisterId<{cName}>
   {{
      public:
         /**
          * @brief Constructor
          */
         {cName}();

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

}} // Backward
}} // Transform
}} // QuICC

#endif // QUICC_TRANSFORM_BACKWARD_{cName.upper()}_HPP
'''
    return content

def source(descr):
    cName = descr['class']
    name = descr['name']
    tag = descr['tag']
    content = f'''/**
 * @file {cName}.cpp
 * @brief Source of the backward transform operator {name}
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/{cName}.hpp"

// Project includes
//

namespace QuICC {{

namespace Transform {{

namespace Backward {{

   std::string {cName}::sTag()
   {{
      return "{tag}";
   }}

   std::string {cName}::sFormatted()
   {{
      return "{name}";
   }}

   {cName}::{cName}()
      : IRegisterId<{cName}>({cName}::sTag(), {cName}::sFormatted())
   {{
   }}

}} // Backward
}} // Transform
}} // QuICC
'''
    return content

names = []
names.append({'class':"D1", 'name':"Backward::D1", 'tag':"Bwd::D1"})
names.append({'class':"D1Laplh", 'name':"Backard::D1Laplh", 'tag':"Bwd::D1Laplh"})
names.append({'class':"D1LaplhZD1R_1D1R1", 'name':"Backard::D1LaplhZD1R_1D1R1", 'tag':"Bwd::D1LaplhZD1R_1D1R1"})
names.append({'class':"D1R1", 'name':"Backard::D1R1", 'tag':"Bwd::D1R1"})
names.append({'class':"D1ZP", 'name':"Backard::D1ZP", 'tag':"Bwd::D1ZP"})
names.append({'class':"D2", 'name':"Backard::D2", 'tag':"Bwd::D2"})
names.append({'class':"D3", 'name':"Backard::D3", 'tag':"Bwd::D3"})
names.append({'class':"DfLaplh", 'name':"Backard::DfLaplh", 'tag':"Bwd::DfLaplh"})
names.append({'class':"DsLaplh", 'name':"Backard::DsLaplh", 'tag':"Bwd::DsLaplh"})
names.append({'class':"Laplh", 'name':"Backard::Laplh", 'tag':"Bwd::Laplh"})
names.append({'class':"LaplhZR_1D1R1", 'name':"Backard::LaplhZR_1D1R1", 'tag':"Bwd::LaplhZR_1D1R1"})
names.append({'class':"P", 'name':"Backard::P", 'tag':"Bwd::P"})
names.append({'class':"P0", 'name':"Backard::P0", 'tag':"Bwd::P0"})
names.append({'class':"R_1", 'name':"Backard::R_1", 'tag':"Bwd::R_1"})
names.append({'class':"R_1D1R1", 'name':"Backard::R_1D1R1", 'tag':"Bwd::R_1D1R1"})
names.append({'class':"R_1LaplhPm", 'name':"Backard::R_1LaplhPm", 'tag':"Bwd::R_1LaplhPm"})
names.append({'class':"R_1Pm", 'name':"Backard::R_1Pm", 'tag':"Bwd::R_1Pm"})
names.append({'class':"R_2", 'name':"Backard::R_2", 'tag':"Bwd::R_2"})
names.append({'class':"SLapl", 'name':"Backard::SLapl", 'tag':"Bwd::SLapl"})
names.append({'class':"SRadLapl", 'name':"Backard::SRadLapl", 'tag':"Bwd::SRadLapl"})
names.append({'class':"Sin_1", 'name':"Backard::Sin_1", 'tag':"Bwd::Sin_1"})
names.append({'class':"Sin_1D1Sin", 'name':"Backard::Sin_1D1Sin", 'tag':"Bwd::Sin_1D1Sin"})
names.append({'class':"Sin_1Dphi", 'name':"Backard::Sin_1Dphi", 'tag':"Bwd::Sin_1Dphi"})
names.append({'class':"Sin_1Laplh", 'name':"Backard::Sin_1Laplh", 'tag':"Bwd::Sin_1Laplh"})
names.append({'class':"Sin_1LaplhDphi", 'name':"Backard::Sin_1LaplhDphi", 'tag':"Bwd::Sin_1LaplhDphi"})

incPath = '../include/QuICC/Transform/Backward/'
srcPath = '../src/Transform/Backward/'
for n in names:
    f = open(incPath+n['class']+'.hpp', 'w')
    f.write(header(n))
    f.close()
    f = open(srcPath+n['class']+'.cpp', 'w')
    f.write(source(n))
    f.close()

