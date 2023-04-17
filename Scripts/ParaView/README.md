# ParaView reader for QuICC visualization data

In order to use the reader, you have two options:

1. Manually add the reader to the plugins using the menu: 
    - `Tools->Manage Plugins`
    - `Load New`
    - Find and select the `PythonQuICCReader.py` module

2. Automatically load the plugin when ParaView is started:
    - Set and export `PV_PLUGIN_PATH` to the directory containing the `PythonQuICCReader.py` module
    - Example in a Linux bash shell:
    ```bash
    export PV_PLUGIN_PATH=/path/to/QuICC/Scripts/ParaView
    ```
    - The export command can be added to your `.bashrc` (or equivalent) to load automatically
