#!/usr/bin/bash

ONAME="visTimeseries_$1_$2.xdmf"

# Create first part of file
cat > ${ONAME} <<EOF
<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
   <Domain>
      <Grid Name="Timeseries" GridType="Collection" CollectionType="Temporal">
EOF

# Loop over range of files
for FID in $(seq -f "%04g" $1 $2);do
   # Build current file name
   FNAME="visState${FID}.hdf5"

   # Get X Resolution
   NX=`h5ls -r ${FNAME} | grep /mesh/grid_x | awk '{print $3}' | sed -e 's/{//' -e 's/}//'`
   # Get Y Resolution
   NY=`h5ls -r ${FNAME} | grep /mesh/grid_y | awk '{print $3}' | sed -e 's/{//' -e 's/}//'`
   # Get Z Resolution
   NZ=`h5ls -r ${FNAME} | grep /mesh/grid_z | awk '{print $3}' | sed -e 's/{//' -e 's/}//'`

   # Get scalar names
   SCALARS=`h5ls -r ${FNAME} | grep "{${NX}, ${NY}, ${NZ}}" | awk '{print $1}' | sed -e 's/\// /' -e 's/\// /' | awk '{print $1}'`
   # Get time
   TIME=`h5dump -y -d '/run/time' ${FNAME} | grep -A 1 'DATA {' | tail -n 1 | awk '{print $1}'`

   # Create grid 
cat >> ${ONAME} <<EOF
         <Grid Name="grid" GridType="Uniform">
            <Time Value="${TIME}" />
            <Topology TopologyType="3DRectMesh" NumberOfElements="${NX} ${NY} ${NZ}"/>
            <Geometry GeometryType="VxVyVz">
               <DataItem Dimensions="${NZ}" NumberType="Float" Precision="8" Format="HDF">
                  ${FNAME}:/mesh/grid_z
               </DataItem>
               <DataItem Dimensions="${NY}" NumberType="Float" Precision="8" Format="HDF">
                  ${FNAME}:/mesh/grid_y
               </DataItem>
               <DataItem Dimensions="${NX}" NumberType="Float" Precision="8" Format="HDF">
                  ${FNAME}:/mesh/grid_x
               </DataItem>
            </Geometry>
EOF
   # Create scalars
   for S in ${SCALARS}; do
cat >> ${ONAME} <<EOF
            <Attribute Name="${S}" AttributeType="Scalar" Center="Node">
               <DataItem Dimensions="${NX} ${NY} ${NZ}" NumberType="Float" Precision="8" Format="HDF">
                  ${FNAME}:/${S}/${S}
               </DataItem>
            </Attribute>
EOF
   done

   # Close grid
cat >> ${ONAME} <<EOF
         </Grid>
EOF

done

# Close file
cat >> ${ONAME} <<EOF
      </Grid>
   </Domain>
</Xdmf>
EOF

