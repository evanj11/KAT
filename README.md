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
- Export plot as PNG and SVG and parameters to .txt files

---

## Installation

### Requirements

- Python 3.9+
- PySide6
- NumPy
- SciPy
- Sympy
- Matplotlib
- Pandas

### Installation Instructions

1. Download the appropriate KAT_{os}.tar.bz2 file
2. Untar the release and the pre-bundled KAT application should be ready for use
3. On MacOS, use the disk image to install KAT directly to Applications

- Anyone desiring to make changes to the either the GUI or analysis scripts can rebundle the app directly using `pyinstaller KAT.spec`

___

## Using KAT

Once `KAT` is installed, operating was designed to be as straightforward as possible:

1. Enter a name for the outputted graphs and parameter information in the `Output Filename` block
2. Upload the data CSV file and the blank CSV file (if present)
3. Select the folder in which you wish to save the outputted data
4. Enter either fluorescence or absorbance from the dropdown menu
   - if absorbance, an additional window will appear where you can input your molar absorptivity 
5. Enter information about the substrate: number of dilutions, dilution factor, and maximum concentration
6. Enter time information, mainly the window with which to compute slopes (default is 10)
   - if you desire, uncheck the `Auto Calculate Linear Range` button to manually select linear range of raw data
7. Finally, select your model it will run!
8. You can now open the graph and kinetic parameters by selecting `Display Graph and Output Values`
___
## Testing

Simulated enzyme data is provided in the `tests` folder, along with information regarding substrate concentrations
Simply download the CSV file, then start KAT, following the steps above
The file named `important_information.txt` in the `tests` folder contains information regarding the substrate concentrations.

