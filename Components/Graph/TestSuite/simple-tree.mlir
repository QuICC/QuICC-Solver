// phys -> mod
// Nr x Nphi x Ntheta
!type_uval = !quiccir.view<3x10x6xf64, "R_DCCSC3D_t">
// L x N x M  ( ..., radial, ...)
!type_umod = !quiccir.view<6x2x6xf64, "C_DCCSC3D_t">

!type_tuval = tensor<3x10x6xf64, "R_DCCSC3D_t">
!type_tumod = tensor<6x2x6xf64, "C_DCCSC3D_t">

// func.func private @simpleTree(%R: tensor<?x?x?xf64>, %Theta: tensor<?x?x?xf64>, %Phi: tensor<?x?x?xf64>) -> (tensor<?x?x?xf64>) {
//   // R
//   %R1 = quiccir.fr.int %R : tensor<?x?x?xf64> -> tensor<?x?x?xf64>
//   %R2 = quiccir.al.int %R1 : tensor<?x?x?xf64> -> tensor<?x?x?xf64>
//   %R3 = quiccir.jw.int %R2 : tensor<?x?x?xf64> -> tensor<?x?x?xf64>
//   // Theta
//   %Th1 = quiccir.fr.int %Theta : tensor<?x?x?xf64> -> tensor<?x?x?xf64>
//   %Th2 = quiccir.al.int %Th1 : tensor<?x?x?xf64> -> tensor<?x?x?xf64>
//   %Th3 = quiccir.jw.int %Th2 : tensor<?x?x?xf64> -> tensor<?x?x?xf64>
//   // Phi
//   %Phi1 = quiccir.fr.int %Phi : tensor<?x?x?xf64> -> tensor<?x?x?xf64>
//   %Phi2 = quiccir.al.int %Phi1 : tensor<?x?x?xf64> -> tensor<?x?x?xf64>
//   %Phi3 = quiccir.jw.int %Phi2 : tensor<?x?x?xf64> -> tensor<?x?x?xf64>
//   // Pol
//   %tmp = quiccir.sub %Th2, %R2 : tensor<?x?x?xf64>, tensor<?x?x?xf64> -> tensor<?x?x?xf64>
//   %Pol = quiccir.add %tmp, %Phi2 : tensor<?x?x?xf64>, tensor<?x?x?xf64> -> tensor<?x?x?xf64>
//   return %Pol : tensor<?x?x?xf64>
// }

func.func private @simpleTree(%R: tensor<?x?x?xf64>, %Theta: tensor<?x?x?xf64>, %Phi: tensor<?x?x?xf64>) -> (tensor<?x?x?xf64>) {
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
  // %0 = tensor.empty() : tensor<?x?x?xf64>
  // %Pol = linalg.add ins(%tmp, %Phi3 : tensor<?x?x?xf64>, tensor<?x?x?xf64>) outs(%0): tensor<?x?x?xf64>) -> tensor<?x?x?xf64>

  return %Pol : tensor<?x?x?xf64>
}

func.func @entry(%thisArr: !llvm.ptr<array<7 x ptr>> {llvm.noalias}, %Polv: !type_umod, %Rv: !type_uval, %Thetav: !type_uval, %Phiv: !type_uval) {
  %R = builtin.unrealized_conversion_cast %Rv : !type_uval to !type_tuval
  %Theta = builtin.unrealized_conversion_cast %Thetav : !type_uval to !type_tuval
  %Phi = builtin.unrealized_conversion_cast %Phiv : !type_uval to !type_tuval
  %RR = tensor.cast %R : !type_tuval to tensor<?x?x?xf64>
  %TTheta = tensor.cast %Theta : !type_tuval to tensor<?x?x?xf64>
  %PPhi = tensor.cast %Phi : !type_tuval to tensor<?x?x?xf64>
  %Polt = call @simpleTree(%RR, %TTheta, %PPhi) : (tensor<?x?x?xf64>, tensor<?x?x?xf64>, tensor<?x?x?xf64>) -> tensor<?x?x?xf64>
  %Pol = tensor.cast %Polt : tensor<?x?x?xf64> to !type_tumod
  /// if this is the only consumer write to existing buffer
  quiccir.materialize %Pol in %Polv : (!type_tumod, !type_umod)
  return
}


// ./bin/quiccir-opt ../examples/simple-tree.mlir --inline --set-quiccir-view-lay --convert-quiccir-to-call --quiccir-view-deallocation --lower-quiccir-alloc --canonicalize --finalize-quiccir-view --convert-func-to-llvm --canonicalize
