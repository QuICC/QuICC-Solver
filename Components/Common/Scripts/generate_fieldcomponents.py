#!/usr/bin/env python

def header(descr, ns):
    cName = descr['class']
    name = descr['name']
    content = f'''/**
 * @file {cName}.hpp
 * @brief {ns} field component {name}
 */

#ifndef QUICC_FIELDCOMPONENTS_{ns.upper()}_{cName.upper()}_HPP
#define QUICC_FIELDCOMPONENTS_{ns.upper()}_{cName.upper()}_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/FieldComponents/{ns}/IRegisterId.hpp"

namespace QuICC {{

namespace FieldComponents {{

namespace {ns} {{

   /**
    * @brief {ns} field component {name}
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

}}
}}
}}

#endif // QUICC_FIELDCOMPONENTS_{ns.upper()}_{cName.upper()}_HPP
'''
    return content

def source(descr, ns):
    cName = descr['class']
    name = descr['name']
    tag = descr['tag']
    content = f'''/**
 * @file {cName}.cpp
 * @brief Source of the {ns} field component {name}
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/FieldComponents/{ns}/{cName}.hpp"

// Project includes
//

namespace QuICC {{

namespace FieldComponents {{

namespace {ns} {{

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

}}
}}
}}
'''
    return content

names = []
names.append({'class':"NotUsed", 'name':"(Not used)", 'tag':"notused"})
names.append({'class':"Phi", 'name':"Phi", 'tag':"phi"})
names.append({'class':"R", 'name':"R", 'tag':"r"})
names.append({'class':"S", 'name':"S", 'tag':"s"})
names.append({'class':"Scalar", 'name':"", 'tag':""})
names.append({'class':"Theta", 'name':"Theta", 'tag':"theta"})
names.append({'class':"X", 'name':"X", 'tag':"x"})
names.append({'class':"Y", 'name':"Y", 'tag':"y"})
names.append({'class':"Z", 'name':"Z", 'tag':"z"})

ns = 'Physical'
incPath = '../include/QuICC/FieldComponents/Physical/'
srcPath = '../src/FieldComponents/Physical/'
for n in names:
    f = open(incPath+n['class']+'.hpp', 'w')
    f.write(header(n, ns))
    f.close()
    f = open(srcPath+n['class']+'.cpp', 'w')
    f.write(source(n, ns))
    f.close()


names = []
names.append({'class':"NotUsed", 'name':"(Not used)", 'tag':"notused"})
names.append({'class':"Phi", 'name':"Phi", 'tag':"phi"})
names.append({'class':"Pol", 'name':"Poloidal", 'tag':"pol"})
names.append({'class':"Q", 'name':"Q", 'tag':"q"})
names.append({'class':"R", 'name':"R", 'tag':"r"})
names.append({'class':"S", 'name':"S", 'tag':"s"})
names.append({'class':"Scalar", 'name':"", 'tag':""})
names.append({'class':"T", 'name':"T", 'tag':"t"})
names.append({'class':"Theta", 'name':"Theta", 'tag':"theta"})
names.append({'class':"Tor", 'name':"Toroidal", 'tag':"tor"})
names.append({'class':"X", 'name':"X", 'tag':"x"})
names.append({'class':"Y", 'name':"Y", 'tag':"y"})
names.append({'class':"Z", 'name':"Z", 'tag':"z"})

ns = 'Spectral'
incPath = '../include/QuICC/FieldComponents/Spectral/'
srcPath = '../src/FieldComponents/Spectral/'
for n in names:
    f = open(incPath+n['class']+'.hpp', 'w')
    f.write(header(n, ns))
    f.close()
    f = open(srcPath+n['class']+'.cpp', 'w')
    f.write(source(n, ns))
    f.close()
