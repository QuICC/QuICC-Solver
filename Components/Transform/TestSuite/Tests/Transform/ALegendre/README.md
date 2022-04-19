Running test with DB data
========================

Metadata filename look like "P_id108_np8_r0_meta.dat"
  - The are for operator "P"
  - Database id is id = 108
  - Mimic data distribution with 8 ranks
  - Mimic the data obtained for rank 0
  - Metadata can be generated using FrameworkLoadSplitterTests

The test is executed as:
./TransformALegendreTests [Poly::P:projector] --id 108 --np 8 --rank 0 --ulp 100
