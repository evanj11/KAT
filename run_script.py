import importlib.util
import traceback
import os
import sys

class StreamRedirector:
    """A simple stream that redirects writes to a callback.
    """

    def __init__(self, callback):
        self.callback = callback

    def write(self, text):
        if text.strip():
            try:
                self.callback(text)
            except Exception:
                pass

    def flush(self):
        pass



def run_internal_script(script_path: str, call_main: bool = True, log_callback=None):
    """
    Dynamically imports and executes a Python script.

    Args:
        script_path (str): Absolute path to the script file.
        call_main (bool): If True, looks for and runs a main() function.
    """
    if not os.path.isfile(script_path):
        if log_callback:
            log_callback(f"[ERROR] Script not found: {script_path}")
        else:
            print(f"[ERROR] Script not found: {script_path}")
        return

    try:
        # Load module from path
        if log_callback:
            sys_stdout_backup = sys.stdout
            sys_stderr_backup = sys.stderr
            sys.stdout = StreamRedirector(log_callback)
            sys.stderr = StreamRedirector(log_callback)

        spec = importlib.util.spec_from_file_location("internal_module", script_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        # Optionally call main()
        if call_main and hasattr(module, 'main') and callable(module.main):
            module.main()

    except Exception as e:
        error_text = traceback.format_exc()
        if log_callback:
            log_callback(f"[ERROR] Failed to run script: {script_path}")
            log_callback(f"[EXCEPTION] {e}")
        else:
            print(error_text)
    finally:
        sys.stdout = sys_stdout_backup
        sys.stderr = sys_stderr_backup

