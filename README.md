# Kinetic Analysis Toolkit (KAT)

`KAT` is a Python-based graphical user interface (GUI) application for analyzing enzyme kinetic data. It supports both classical models( **Michaelis-Menten** and **Hill** 
models) and complex models (**Monod-Wyman-Changeux** and **Koshland-Nemethy-Filmer** models), outputting kinetic parameters (*e.g.* *K<sub>M</sub>*, *V<sub>max</sub>*, and *k<sub>cat</sub>* for Michaelis-Menten) along 
with visual plots of your data. KAT is designed for experimental scientists who want a quick, user-friendly tool to interpret enzyme kinetics from absorbance or fluorescence 
data stored in CSV files.

---

## Features

- GUI-based interface (built with PySide6)
- Input CSV files containing absorbance or fluorescence data
- Can accept both replica and mutant data, treating each differently depending on the application
- Supports **Michaelis-Menten**, **Hill**, **Monod-Wyman-Changeux**, and **Koshland-Nemethy-Filmer** kinetic models
- Outputs:
  - Kinetic parameters based on the model chosen including, *V<sub>max</sub>*, *K<sub>M</sub>*, *n*, and many others for complex models
- Automatically fits nonlinear models using SciPy and Sympy
- Displays graph chosen model
- Export plot as both PNG and SVG and parameters to .txt files

---

## Installation

### Requirements

- Python 3.12+
- PySide6
- Matplotlib
- NumPy
- Pandas
- Scikit-learn
- SciPy
- Sympy

### Installation Instructions

1. Download the appropriate KAT_{os}.tar.bz2 file from `Releases`
   - [KAT for MacOS](https://github.com/evanj11/KAT/releases/download/v1.0.0/KAT.dmg)$^\dagger$
   - [KAT for Windows (under development)](https://github.com/evanj11/KAT/releases/download/v1.0.0/Kinetic.Analysis.Toolkit-0.0.1.msi)
3. Untar the release and the pre-bundled KAT application should be ready for use
4. On MacOS, use the disk image to install KAT directly to Applications

- Anyone desiring to make changes to the either the GUI or analysis scripts can rebundle the app directly using `pyinstaller KAT.spec`
- Those wishing to run KAT outside of a dedicated app can run `python main.py` in the master directory to access the GUI<br />

$^\dagger$ currently only supports Apple Silicon chips (M1 and later)<br />
  - For Windows and Linux users:
     
      ```
      git clone https://github.com/evanj11/KAT
      cd KAT
      pip install -r requirements.txt
      python main.py
      ```
___
    
## Using KAT

### Running through the GUI

Once `KAT` is installed, operating was designed to be as straightforward as possible:

1. Enter a name for the outputted graphs and parameter information in the `Output Filename` block
2. Upload the data CSV file and the blank CSV file (if present)
   - Make sure data is formatted like the following:
  
     | Time Data | Temperature | Max Substrate | Next Subsubstrate | ... |
     | :---: | :---: | :---: | :---: | :---: |
     | 00:00:00 |  29C  | raw data   | raw data | ... |
     | 00:01:00 |  29C  | raw data   | raw data | ... |
     | ... | ... | ... | ... | ... |

4. Select the folder in which you wish to save the outputted data
5. Enter either fluorescence or absorbance from the dropdown menu
   - if absorbance, an additional window will appear where you can input your molar absorptivity 
6. Enter information about the substrate: number of dilutions, dilution factor, and maximum concentration
7. Enter time information, mainly the window with which to compute slopes (default is 10)
   - if you desire, uncheck the `Auto Calculate Linear Range` button to manually select linear range of raw data
8. Finally, select the model you want to run!
9. You can now open the graph and kinetic parameters by selecting `Display Graph and Output Values`

### Running through the terminal

If you do not want to use the GUI or already have the velocity information, `KAT` can be run through the command line:

    python scripts/kat_cmd.py

Running 
  ```
  python scripts/kat_cmd.py --help
  ```
will open the help menu where all command-line arguments can be listed. The required arguments are:
  1. -s/--substrate = either a list or file containing substrate concentrations or the number of substrate concentrations, the dilution factor, and the maximum substrate concentration
     - if the later is supplied, use the --calc-sub-dil flag to automatically generate the substrate concentrations from the dilution information provided
  2. -v/--velocity = either a list or file containing the velocities at each substrate concentration
  3. -m/--model = the chosen kinetic model to run chosen from `mm` (*M*ichaelis-*M*enten), `hill` (Hill Kinetics), `mwc` (*M*onod-*W*yman-*C*hangeux), or `knf` (*K*oshland-*N*emethy-*F*ilmer)
     - if `knf` model chosen, the --sites flag must be used to specify the number of catalytic pockets in the enzyme
  5. -o/--output = the output filename 
___
### Example Usage of Command Line KAT

If providing a list of substrate and velocity values and computing Michaelis-Menten Kinetics:
  ```
  python scripts/kat_cmd.py -s 10 20 35 50 65 80 100 -v 0.0566 0.085 0.1285 0.1478 0.1662 0.1774 0.1799 -m mm
  ```

If providing files for substrate and velocity data:
  ```
  python scripts/kat_cmd.py -s substrate.txt -v velocities.txt -m knf --calc-sub-dil --sites 2
  ```
and substrate.txt looks like:
  ```
  7     #number of dilutions
  2     #dilution factor
  100   #maximum substrate value
  ```

## Testing

Simulated enzyme data is provided in the `tests` folder, along with information regarding substrate concentrations.
Simply download the CSV file, then start KAT, following the steps above. 
The file named `important_information.txt` in the `tests` folder contains information regarding the substrate concentrations.<br />


