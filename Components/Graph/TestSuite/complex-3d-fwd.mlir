// phys -> mod
// Nr x Nphi x Ntheta
// !type_uval = !quiccir.view<3x10x6xf64, "R_DCCSC3D_t">
// // L x N x M  ( ..., radial, ...)
// !type_umod = !quiccir.view<6x2x6xf64, "C_DCCSC3D_t">

func.func @entry(%R: tensor<?x?x?xf64>, %Theta: tensor<?x?x?xf64>, %Phi: tensor<?x?x?xf64>) -> (tensor<?x?x?xf64>) {
  // R
  %R1 = quiccir.fr.int %R : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 0 :i64}
  %R1T = quiccir.transpose %R1 permutation = [2, 0, 1] : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 1 :i64}
  %R2 = quiccir.al.int %R1T : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 2 :i64}
  %R2T = quiccir.transpose %R2 permutation = [2, 0, 1] : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 3 :i64}
  %R3 = quiccir.jw.int %R2T : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 4 :i64}
  // Theta
  %Th1 = quiccir.fr.int %Theta : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 0 :i64}
  %Th1T = quiccir.transpose %Th1 permutation = [2, 0, 1] : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 1 :i64}
  %Th2 = quiccir.al.int %Th1T : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 2 :i64}
  %Th2T = quiccir.transpose %Th2 permutation = [2, 0, 1] : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 3 :i64}
  %Th3 = quiccir.jw.int %Th2T : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 4 :i64}
  // Phi
  %Phi1 = quiccir.fr.int %Phi : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 0 :i64}
  %Phi1T = quiccir.transpose %Phi1 permutation = [2, 0, 1] : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 1 :i64}
  %Phi2 = quiccir.al.int %Phi1T : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 2 :i64}
  %Phi2T = quiccir.transpose %Phi2 permutation = [2, 0, 1] : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 3 :i64}
  %Phi3 = quiccir.jw.int %Phi2T : tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 4 :i64}
  // Pol
  %tmp = quiccir.sub %Th3, %R3 : tensor<?x?x?xf64>, tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 5 :i64}
  %Pol = quiccir.add %tmp, %Phi3 : tensor<?x?x?xf64>, tensor<?x?x?xf64> -> tensor<?x?x?xf64> attributes{implptr = 6 :i64}
  return %Pol : tensor<?x?x?xf64>
}

// /home/gcastigl/codes/quiccir/build-17/bin/quiccir-opt ./complex-3d-fwd.mlir --inline --quiccir-view-wrapper='dim-rets=7,3,6 dim-args=6,10,10 lay-args=R_DCCSC3D_t lay-rets=C_DCCSC3D_t' --inline --inline --set-quiccir-dims='phys=6,10,10 mods=3,6,7' --set-quiccir-view-lay='lay-ppp2mpp=R_DCCSC3D_t,C_DCCSC3D_t lay-pmp2mmp=C_DCCSC3D_t,C_S1CLCSC3D_t lay-pmm2mmm=C_DCCSC3D_t,C_DCCSC3D_t' --canonicalize --convert-quiccir-to-call --quiccir-view-deallocation --lower-quiccir-alloc --canonicalize --finalize-quiccir-view --cse --convert-func-to-llvm --canonicalize
