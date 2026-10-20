# Debugging the MPI version of QuICC with GDB

## hook the running process to gdb

In the run folder:

        QUICC_GDB_HOOK=0 mpirun -n 4 /QuICC/QuICC-Solver-merge/build_NL_anelastic_mpi_debug/Models/AnelasticSphereRTC/Executables/AnelasticSphereRTCExplicitHydroModel

The output looks like:

        PID 84255 on d429e112fba6 with rank 0 ready to attach

Open another terminal in the same container and start gdb

        gdb -p 84255


## breakout lines (example)

                (gdb) break IPhysicalKernel.cpp:86

                (gdb) bt

Find the gdbHook line that looks like

                #3  0x0000aaaacd581ff8 in QuICC::Environment::Mpi::gdbHook (this=0xaaaace780858 <QuICC::QuICCEnv()::env>) at /QuICC/QuICC-Solver-merge/Components/Environment/Environment/Mpi.cpp:192

Then:

                (gdb) frame 3

To see the local variables:

                (gdb) info locals

Inspect loop

                (gdb) list Mpi.cpp:180,200

Change wait condition and continue:

                (gdb) set var wait = false
                (gdb) continue

Now we are at the breakpoint

                (gdb) set print pretty on
                
Then we can print some variables

                (gdb) p this->mVectors.size()
                (gdb) p name

Print the call stack

                (gdb) bt

To move to the first frame of the stack:

                (gdb) frame 1
