import sys
import subprocess

def main():
    subprocess.run([sys.executable, "kat/toolkit/kat_gui-ctk.py"], check=True)
