# Kinetic Analysis Toolkit (KAT)

`KAT` is a Python-based graphical user interface (GUI) application for analyzing enzyme kinetic data. It supports both classical models( **Michaelis-Menten** and **Hill** 
models) and complex models (**Monod-Wyman-Changeux** and **Koshland-Nemethy-Filmer** models), outputting kinetic parameters (Km, Vmax, kcat for Michaelis-Menten) along 
with visual plots of your data. KAT is designed for experimental scientists who want a quick, user-friendly tool to interpret enzyme kinetics from absorbance or fluorescence 
data stored in CSV files.

---

## Features

- GUI-based interface (built with PySide6)
- Input CSV files containing absorbance or fluorescence data
- Can accept both replica and mutant data, treating each differently depending on the application
- Supports **Michaelis-Menten**, **Hill**, **Monod-Wyman-Changeux**, and **Koshland-Nemethy-Filmer** kinetic models
- Outputs:
  - Km (Michaelis constant)
  - Vmax (maximum velocity)
  - kcat (turnover number)
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

Download the appropriate KAT_{os}.tar.bz2 file
Untar the release and the pre-bundled KAT application should be ready for use
On MacOS, use the disk image to install KAT directly to Applications

- Anyone desiring to make changes to the either the GUI or analysis scripts can rebundle the app directly using `pyinstaller KAT.spec`

___

## Using KAT



___
## Testing

Simulated enzyme data is provided in the `tests` folder, along with information regarding substrate concentrations
Simply download the CSV file, then start KAT

