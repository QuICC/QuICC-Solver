# Debugging in VScode

## Connect VScode to running container

We assume we are running in a container.

- Start VScode

- Install C++ extension

- Istall Dev Containers

- Start the docker container in a terminal window. For example, this containeris for the nonlinear models and will run in the `~/QuICC/QuICC_mnt` folder :

        docker run --rm -it \
        -v ~/QuICC/QuICC_mnt:/QuICC \
        --mount type=bind,src="/run/host-services/ssh-auth.sock",target=/run/host-services/ssh-auth.sock \
        -e SSH_AUTH_SOCK="/run/host-services/ssh-auth.sock" \
        quicc/buildbase:latest


- Attach VScode to running container. From the page: https://code.visualstudio.com/docs/devcontainers/attach-container

    Command palette -> F1 -> Dev Containers: Attach to Running Container -> select container (<funny name> quicc/buildbase:latest )

    Internet connection might be needed for this first step. But the session should keep working if internet is lost, **provided the docker session to which vscode is now attached is not closed.**

- Re-install the following extension

    C/C++

    C/C++ Extension Pack

    Code Runner (usually recommended, but I don't think it's needed)

- Open project folder. In our case `/QuICC/QuICC/`


## Connect to remote host (eg LUMI)

Follow instructions on https://code.visualstudio.com/docs/remote/ssh 

Specifically:
- Install Remote-SSH extension
- in the command palette of VSCode:

        >Remote-SSH: Connect to Host
- Add new host if it's the first time: maffeist@lumi.cscs.fi
- Choose the `.ssh/config` file to login. *things need to be set up properly for a ssh login on this host*
- Enter password


## Debug


- If needed, install the gdb in the running container (from the terminal window):

        apt-get update
        apt-get install gdb

- Install the VScode extension `GDB Debugger - Beyond`

- Compile QuICC in debug mode. For example:
    
        mkdir build_NL_anelastic_debug

        cd build_NL_anelastic_debug

        # to compile the nonlinear version of the model

        cmake .. -DQUICC_MODEL=AnelasticShellRTC -DCMAKE_BUILD_TYPE=Debug

        make -j 6

        # to compile the linear stability version

        cmake .. -DQUICC_MODEL=AnelasticShellRTC -DCMAKE_BUILD_TYPE=Debug -DPETSC_DIR=/opt/view -DSLEPC_DIR=/opt/view

        make -j 6 AnelasticShellRTCImplicitStability
        
        # Create a folder where the debug is to be run
        # needs to be a place visible in the container
        # example

        mkdir /QuICC/AnelasticShellRTC/debug
        cd /QuICC/AnelasticShellRTC/debug

        # assuming we want to debug the Model executable:
        # create/copy a parameters.cfg and state_initial.hdf5

- Start the debugger as shown on https://marketplace.visualstudio.com/items?itemName=coolchyni.beyond-debug :

    - Click the run and debug button
    - create a launch.json file

- Copy the following in the launch.json file (which is save in .vscode):

        {
            // Use IntelliSense to learn about possible attributes.
            // Hover to view descriptions of existing attributes.
            // For more information, visit: https://go.microsoft.com/fwlink/?linkid=830387
            "version": "0.2.0",
            "configurations": [

                {
                    "type": "by-gdb",
                    "request": "launch",
                    "name": "Launch(gdb)",
                    //"program": "${fileBasenameNoExtension}",
                    "program": "/QuICC/QuICC/build_NL_anelastic_debug/Models/AnelasticShellRTC/Executables/AnelasticShellRTCExplicitModel",
                    "cwd": "/QuICC/AnelasticShellRTC/debug/"
                }
            ]
        }

- click the launch green button in the top left corner

- useful instruction for debugging in VScode can be found, for example:
    - https://code.visualstudio.com/docs/editor/debugging


## Debug tricks


### Kill hanging run

Sometimes, you press the "run" button and VSCode hangs there. If there is the process bar in the "run and debug" panel, you can kill that process with Shift-F5 (Shift-fn-F5 on a Mac)


### Print variable values

In the debug console, once a breakpoint is hit, the command `print` can be used to print variables values. 

Clearly the breakpoint needs to be in a sensible place.

For example to print an int (example with a breakpoint in `MomentumKernel.cpp`):

    print this->mInertia

    $1 = 1

### Print arrays

#### Example 1: mRadius

Example with a breakpoint in `MomentumKernel.cpp`:

    print this->mRadius

    $17 = {<Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >> = {<Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >> = {<Eigen::DenseBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >> = {<Eigen::DenseCoeffsBase<Eigen::Matrix<double, -1, 1, 0, -1, 1>, 3>> = {<Eigen::DenseCoeffsBase<Eigen::Matrix<double, -1, 1, 0, -1, 1>, 1>> = {<Eigen::DenseCoeffsBase<Eigen::Matrix<double, -1, 1, 0, -1, 1>, 0>> = {<Eigen::EigenBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >> = {<No data fields>}, <No data fields>}, <No data fields>}, <No data fields>}, <No data fields>}, <No data fields>}, m_storage = {m_data = 0xaaaaad160540, m_rows = 63}}, <No data fields>}

which is basically the information that one can desume from the Variables->Locals panel. However, we get also the size of the array (64).

    print this->mRadius.data()[0]@10
    
    $19 = {1.5383061294619464, 1.5370634370521283, 1.5345811417616244, 1.5308654149676424, 1.5259254945524501, 1.5197736619365443, 1.5124252115451042, 1.5038984127836403, 1.4942144646173747, 1.4833974428672727}

To print the first 10 values of `mRadius`. And in this case.

    print this->mRadius.data()[0]@63

prints the whole array.


#### Example 2: rNLComp

The above depends on the available methods in a specific variable.

`profile` gives a 1D array

    print rNLComp.profile(1,1).m_data[0]@10

    $5 = {-5.7529597132538956e-22, -3.4987299068798639e-22, -1.8358961327007484e-22, -8.9977341483363664e-23, -6.7771018307778425e-23, -1.0163276570788058e-22, -1.6797163942900139e-22, -2.4280020882199813e-22, -3.0791645186783028e-22, -3.5333261459505452e-22}

`slice` gives a 2D matrix

    print rNLComp.slice(1).data()[0]@10
    
    $10 = {-7.0422654329261449e-22, -6.2139287690149649e-22, -5.607314587353385e-22, -5.2212513267841266e-22, -5.0383597814835437e-22, -5.0281346547060178e-22, -5.1514657175950327e-22, -5.365816772684157e-22, -5.6301657379493283e-22, -5.9089395740012698e-22}

#### Example 3: print array sizes:

for arrays:
    
    print rho.size()

for matrices

    print tmpSquare.cols()
    print tmpSquare.rows()

