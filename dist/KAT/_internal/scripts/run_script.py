import importlib.util
import traceback
import os
import sys

def run_internal_script(script_path: str, call_main: bool = True):
    """
    Dynamically imports and executes a Python script.

    Args:
        script_path (str): Absolute path to the script file.
        call_main (bool): If True, looks for and runs a main() function.
    """
    if not os.path.isfile(script_path):
        print(f"[ERROR] Script not found: {script_path}")
        return

    try:
        # Load module from path
        spec = importlib.util.spec_from_file_location("internal_module", script_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        # Optionally call main()
        if call_main and hasattr(module, 'main') and callable(module.main):
            module.main()

    except Exception:
        print(f"[ERROR] Failed to run script: {script_path}")
        traceback.print_exc()

