#!/usr/bin/bash

NX=`h5ls -r $1 | grep /mesh/grid_x | awk '{print $3}' | sed -e 's/{//' -e 's/}//'`
NY=`h5ls -r $1 | grep /mesh/grid_y | awk '{print $3}' | sed -e 's/{//' -e 's/}//'`
NZ=`h5ls -r $1 | grep /mesh/grid_z | awk '{print $3}' | sed -e 's/{//' -e 's/}//'`

FID=`basename -s .hdf5 $1 | sed -e 's/visState//'`
SCALARS=`h5ls -r $1 | grep "{${NX}, ${NY}, ${NZ}}" | awk '{print $1}' | sed -e 's/\// /' -e 's/\// /' | awk '{print $1}'`
TIME=`h5dump -y -d '/run/time' $1 | grep -A 1 'DATA {' | tail -n 1 | awk '{print $1}'`

cat > vis${FID}.xdmf <<EOF
<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
   <Domain>
      <Grid Name="grid" GridType="Uniform">
         <Topology TopologyType="3DRectMesh" NumberOfElements="${NX} ${NY} ${NZ}"/>
         <Geometry GeometryType="VxVyVz">
            <DataItem Dimensions="${NZ}" NumberType="Float" Precision="8" Format="HDF">
               visState${FID}.hdf5:/mesh/grid_z
            </DataItem>
            <DataItem Dimensions="${NY}" NumberType="Float" Precision="8" Format="HDF">
               visState${FID}.hdf5:/mesh/grid_y
            </DataItem>
            <DataItem Dimensions="${NX}" NumberType="Float" Precision="8" Format="HDF">
               visState${FID}.hdf5:/mesh/grid_x
            </DataItem>
         </Geometry>
EOF
for S in ${SCALARS}; do
cat >> vis${FID}.xdmf << EOF
         <Attribute Name="${S}" AttributeType="Scalar" Center="Node">
            <DataItem Dimensions="${NX} ${NY} ${NZ}" NumberType="Float" Precision="8" Format="HDF">
               visState${FID}.hdf5:/${S}/${S}
            </DataItem>
         </Attribute>
EOF
done
cat >> vis${FID}.xdmf <<EOF
         <Time Value="${TIME}" />
      </Grid>
   </Domain>
</Xdmf>
EOF

