"""Lets `submit --local` start papermill from a conda environment such as nanoHUB's TorchSim.

simtool runs `submit --local <env>/bin/papermill ...`. submit starts it without the environment's library path,
so its compiled modules load the system libstdc++ (Debian 10), which is too old for the environment's libzmq:
    ImportError: ... libstdc++.so.6: version `GLIBCXX_3.4.29' not found (required by ... libzmq.so.5)
and papermill never runs. The Jupyter kernel itself is fine because it is started with that path.

use_env_libraries() writes a small wrapper next to the caller that puts <env>/lib first on LD_LIBRARY_PATH and
then runs the environment's papermill, and tells simtool to use it. Off nanoHUB, or outside a conda
environment, it does nothing.
"""
import sys
from pathlib import Path


def use_env_libraries(where="."):
    import simtool.run
    env = Path(sys.prefix)
    papermill, lib = env / "bin" / "papermill", env / "lib"
    if not simtool.run.submitAvailable or not papermill.exists() or not (lib / "libstdc++.so.6").exists():
        return None
    wrapper = Path(where, ".papermill_env").resolve()        # absolute: submit runs in the run folder
    wrapper.write_text("#!/bin/sh\n"
                       f'LD_LIBRARY_PATH="{lib}${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}"\n'
                       "export LD_LIBRARY_PATH\n"
                       f'exec "{papermill}" "$@"\n')
    wrapper.chmod(0o755)
    simtool.run.papermillCLI = str(wrapper)                  # read when each Run starts
    return wrapper
