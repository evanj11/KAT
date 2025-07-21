#!/usr/bin/env python
# coding: utf-8

# In[5]:


import sys
import os
import subprocess
import runpy
import traceback
from PySide6.QtWidgets import (
    QApplication, QWidget, QLabel, QLineEdit, QPushButton,
    QTabWidget, QFileDialog, QVBoxLayout, QHBoxLayout, QTextEdit,
    QGridLayout, QMainWindow, QTreeWidget, QTreeWidgetItem, QComboBox, 
    QGroupBox, QSizePolicy, QGraphicsOpacityEffect, QMessageBox, QCheckBox, QStyle
)
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtGui import QPixmap, QFont, QPainter, QPalette, QColor, QImage, QDoubleValidator, QMovie
from PySide6.QtCore import Qt, QTimer, QSize, QThread, QObject, Signal

import tempfile
import shutil

def get_movie_path(relative_path):
    if getattr(sys, 'frozen', False):
        base_path = sys._MEIPASS
    else:
        base_path = os.path.abspath(os.path.dirname(__file__))
    resource_path = os.path.join(base_path, relative_path)

    # Extract to a temp file because QMovie needs a file on disk
    temp_dir = tempfile.gettempdir()
    temp_path = os.path.join(temp_dir, os.path.basename(relative_path))
    shutil.copyfile(resource_path, temp_path)
    return temp_path


def get_resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller/briefcase. """
    if getattr(sys, 'frozen', False):
	    base_path = sys._MEIPASS
    else:
	    base_path = os.path.abspath(os.path.dirname(__file__))
    return os.path.join(base_path, relative_path)

SCRIPTS_DIR = get_resource_path('scripts')

if SCRIPTS_DIR not in sys.path:
    sys.path.append(SCRIPTS_DIR)


def safe_open_file_dialog(parent, caption="Open File", file_filter="CSV Files (*.csv);;All Files (*)"):
    try:
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(
            parent, caption, "", file_filter, options=options
        )
        if file_path:
            return file_path
        else:
            # Ask user for manual path input if they cancel or fail to select
            fallback, ok = QInputDialog.getText(
                parent, "Manual File Entry",
                "Enter full path to CSV file:"
            )
            if ok and os.path.isfile(fallback):
                return fallback
            elif ok:
                QMessageBox.warning(parent, "Invalid Path", "File not found at the specified location.")
    except Exception as e:
        print("File dialog error:", e)
        traceback.print_exc()
        QMessageBox.critical(parent, "Error", "Could not open file dialog due to an internal error.")
    return None


class ScriptWorker(QObject):
    finished = Signal()

    def __init__(self, script_path):
        super().__init__()
        self.script_path = script_path

    def run(self):
        from run_script import run_internal_script  # Import here to avoid freezing
        run_internal_script(self.script_path, call_main=True)
        self.finished.emit()


class GeneralHelpWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Help")
        self.setGeometry(200, 200, 400, 300)
        layout = QVBoxLayout()
        help_text = QLabel(
            "Welcome to KAT, the Kinetic Analysis Toolkit!\n"
            "To generate a kinetic curve following either Hill or Michelis-Menten kinetics, follow these steps:\n"
            "1. Create a .csv file with values for fluorescence or absorbance that follows [Time] [Temp] [Data].\n"
            "2. Give the output graph file a name.\n"
            "3. Specify a path to the data in a .csv format.\n"
            "4. Specify the kind of data you have collected (absorbance vs. fluorescence).\n"
            "5. Input the information about the substrate concentrations used (must be serially diluted).\n"
            "6. Input the information about the time-course to sample.\n"
            "        Warning: Time step must evenly divide the difference between Time Min and Time Max \n"
            "        (the quotient being the minutes between slope calculations).\n"
            "7. Select the columns to read data from.\n"
            "8. Click either *Graph Averages* or *Graph Best-Fit* for Hill or Michelis-Menten Kinetics.\n"
            "9. Click *Display Graph* to open the graph output."
        )
        help_text.setFont(QFont("Times New Roman", 14))
        layout.addWidget(help_text)
        self.setLayout(layout)

class TimeHelpWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Help with Inputting Time Information")
        self.setGeometry(200, 200, 400, 300)
        layout = QVBoxLayout()
        help_text = QLabel(
            "To auto-compute the linear range of the velocities, leave the Auto Compute Linear Range button clicked.\n"
            "This will compare the percent difference between two subsequent slopes and begin the linear range when\n"
            "the percent difference between two slopes is less than 5% (on average); however, if you would like to\n"
            "input your own start and stop time for the linear range, you can unclick the Auto Compute Linear Range button.\n"
            "Then, you can input the minimum and maximum time to compute velocities and the number of velocities to compute\n"
            "(i.e. the Time Step). The Time Step must be able to divide the difference in Time Max and Time Min.\n"
            "Regardless of using the Auto Compute Linear Range function, you must input the window to compute velocities\n"
            "(i.e. the number of data points between the computing of a slope)."
        )
        help_text.setFont(QFont("Times New Roman", 14))
        layout.addWidget(help_text)
        self.setLayout(layout)

class DisplayComplexData(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Complex Data")
        self.setGeometry(400, 100, 200, 200)
        layout = QVBoxLayout()
        work_dir = os.environ.get('WORKING_DIR')
    
        with open(os.path.join(work_dir, 'name_data.txt'), 'r') as f:
            name = [line.strip() for line in f.readlines()]

        with open(os.path.join(work_dir, f'{name[0]}_complex_data.txt'), 'r') as f:
            complex_data = [line.strip() for line in f.readlines()]
        if 'Poor confidence' in complex_data:
            QMessageBox.warning(self, 'Warning', 'Complex model resulted in low confidence results, proceed with caution')

        if 'Monod-Wyman-Changeux' in complex_data:
            layout.addWidget(QLabel("Monod-Wyman-Changeux Results:\n"
                                    f"    Velocity Tensed State = {complex_data[2]}\n"
                                    f"    Velocity Relaxed State = {complex_data[3]}\n"
                                    f"    K in Tensed State = {complex_data[4]}\n"
                                    f"    K in Relaxed State = {complex_data[5]}\n"
                                    f"    Allosteric Constant = {complex_data[6]}\n"
                                    f"    Number of Allosteric Sites = {complex_data[7]}\n"
                                    ))
            if 'Using Cross-Validation values' in complex_data:
                layout.addWidget(QLabel("\nCross-Validation values are displayed"))
            else:
                layout.addWidget(QLabel("\nStandard Best-Fit values are displayed"))
        elif 'Koshland-Nemethy-Filmer' in complex_data:
            layout.addWidget(QLabel(f"Koshland-Nemethy-Filmer Results:\n"
                                    f"    Vmax = {complex_data[2]}\n"
                                    f"    Kd = {complex_data[3]}\n"
                                    f"    Basal Activity = {complex_data[4]}\n"
                                    f"    Gamma = {complex_data[5]}\n"
                                    ))
            if 'Using Cross-Validation values' in complex_data:
                layout.addWidget(QLabel("\nCross-Validation values are displayed"))
            else:
                layout.addWidget(QLabel("\nStandard Best-Fit values are displayed"))
        self.setLayout(layout)



class DisplayMutData(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Mutant Data")
        self.setGeometry(400, 100, 200, 200)
        layout = QVBoxLayout()
        work_dir = os.environ.get('WORKING_DIR')

        with open(os.path.join(work_dir, 'name_data.txt'), 'r') as f:
            name = [line.strip() for line in f.readlines()]

        mut_data = []
        with open(os.path.join(work_dir, f'{name[0]}_mutant_data.txt'), 'r') as f:
            mut_data_temp = [line.strip() for line in f.readlines()]
            for mut in mut_data_temp:
                mut_split = mut.split(",")
                mut_data.append(mut_split)
        if os.path.exists(os.path.join(work_dir, f'{name[0]}_complex_data.txt')):
            with open(os.path.join(work_dir, f'{name[0]}_complex_data.txt'), 'r') as f:
                complex_data = [line.strip() for line in f.readlines()]
        else:
            complex_data = []
        if 'Poor confidence' in complex_data:
            QMessageBox.warning(self, 'Warning', 'Complex model resulted in low confidence results, proceed with caution')

        for i in range(len(mut_data)):
            mutant = mut_data[i]
            if 'Monod-Wyman-Changeux' in complex_data:
                if i == 0:
                    layout.addWidget(QLabel("Wild-Type Results:\n" 
                                            "Monod-Wyman-Changeux Results:\n"
                                            f"    Velocity Tensed State = {mutant[0]}\n"
                                            f"    Velocity Relaxed State = {mutant[1]}\n"
                                            f"    K in Tensed State = {mutant[2]}\n"
                                            f"    K in Relaxed State = {mutant[3]}\n"
                                            f"    Allosteric Constant = {mutant[4]}\n"
                                            f"    Number of Allosteric Sites = {mutant[5]}\n"
                                            ))
                else:
                    layout.addWidget(QLabel(f"Mutant #{i} Results:\n"
                                            "Monod-Wyman-Changeux Results:\n"
                                            f"    Velocity Tensed State = {mutant[0]}\n"
                                            f"    Velocity Relaxed State = {mutant[1]}\n"
                                            f"    K in Tensed State = {mutant[2]}\n"
                                            f"    K in Relaxed State = {mutant[3]}\n"
                                            f"    Allosteric Constant = {mutant[4]}\n"
                                            f"    Number of Allosteric Sites = {mutant[5]}\n"
                                            ))
                if 'Using Cross-Validation values' in complex_data:
                    layout.addWidget(QLabel("Cross-Validation values are displayed\n"))
                else:
                    layout.addWidget(QLabel("Standard Best-Fit values are displayed\n"))

            elif 'Koshland-Nemethy-Filmer' in complex_data:
                if i == 0:
                    layout.addWidget(QLabel("Wild-Type Results:\n"
                                        "Koshland-Nemethy-Filmer Results:\n"
                                        f"    Vmax = {mutant[0]}\n"
                                        f"    Kd = {mutant[1]}\n"
                                        f"    Basal Activity = {mutant[2]}\n"
                                        f"    Gamma = {mutant[3]}\n"
                                        ))
                else:
                    layout.addWidget(QLabel(f"Mutant #{i} Results:\n"
                                        f"Koshland-Nemethy-Filmer Results:\n"
                                        f"    Vmax = {mutant[0]}\n"
                                        f"    Kd = {mutant[1]}\n"
                                        f"    Basal Activity = {mutant[2]}\n"
                                        f"    Gamma = {mutant[3]}\n"
                                        ))
                if 'Using Cross-Validation values' in complex_data:
                    layout.addWidget(QLabel("Cross-Validation values are displayed\n"))
                else:
                    layout.addWidget(QLabel("Standard Best-Fit values are displayed\n"))

            else:
                if i == 0:
                    layout.addWidget(QLabel("Wild-Type Results:\n"
                                            f"    Hill Coefficient = {mutant[0]}\n"
                                            f"    Vmax = {mutant[1]}\n"
                                            f"    Km = {mutant[2]}\n"
                                            ))
                else:
                    layout.addWidget(QLabel(f"Mutant #{i} Results:\n"
                                            f"    Hill Coefficient = {mutant[0]}\n"
                                            f"    Vmax = {mutant[1]}\n"
                                            f"    Km = {mutant[2]}\n"
                                            ))
        self.setLayout(layout)

class BatchImportCSV(QWidget):
    folder_selected = Signal(str)
    done = Signal()
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Batch Import CSV Files for Replica Analysis")
        self.setGeometry(100, 200, 200, 400)
        layout = QGridLayout()
        self.batch_import_label = QLabel("Please select how many CSV files to import\n"
                                         "Note: for >5 replicas, please repeat this process")
        layout.addWidget(self.batch_import_label, 0, 0)
        self.num_rep_combo = QComboBox()
        self.num_rep_combo.setMaximumWidth(150)
        self.num_rep_combo.addItems(["Select", "2", "3", "4", "5"])
        layout.addWidget(QLabel("Number of Replicas to Analyze:"), 1, 0)
        layout.addWidget(self.num_rep_combo, 1, 1)
        self.mut_check = QCheckBox("Are Mutants Present?")
        self.mut_check.setChecked(False)  # Optional: set initial state
        layout.addWidget(self.mut_check, 1, 2) 
        self.mut_check.stateChanged.connect(self.mut_toggle)        
        self.num_rep_combo.currentTextChanged.connect(self.toggle_data_rep_input)
       
        self.rep1_path_input = QLineEdit()
        self.rep1_path_input.setMaximumWidth(225)
        self.rep1_path_input.hide()
        self.rep1_path_button = QPushButton("Data CSV")
        self.rep1_path_button.setMaximumWidth(100)
        self.rep1_path_button.clicked.connect(self.select_rep1_input_file)
        self.rep1_path_button.hide()
        self.rep1_label = QLabel("Replica #1 Data CSV File Path:")
        self.rep1_label.hide()
        self.wt_label = QLabel("Wild-Type Data CSV File Path:")
        self.wt_label.hide()
        layout.addWidget(self.wt_label, 2, 0)
        layout.addWidget(self.rep1_label, 2, 0)
        layout.addWidget(self.rep1_path_input, 2, 1)
        layout.addWidget(self.rep1_path_button, 2, 2)
        
        self.rep2_path_input = QLineEdit()
        self.rep2_path_input.setMaximumWidth(225)
        self.rep2_path_input.hide()
        self.rep2_path_button = QPushButton("Data CSV")
        self.rep2_path_button.setMaximumWidth(100)
        self.rep2_path_button.clicked.connect(self.select_rep2_input_file)
        self.rep2_path_button.hide()
        self.rep2_label = QLabel("Replica #2 Data CSV File Path:")
        self.rep2_label.hide()
        self.mut2_label = QLabel("Mutant #1 Data CSV File Path:")
        self.mut2_label.hide()
        layout.addWidget(self.mut2_label, 3, 0)
        layout.addWidget(self.rep2_label, 3, 0)
        layout.addWidget(self.rep2_path_input, 3, 1)
        layout.addWidget(self.rep2_path_button, 3, 2)

        self.rep3_path_input = QLineEdit()
        self.rep3_path_input.setMaximumWidth(225)
        self.rep3_path_input.hide()
        self.rep3_path_button = QPushButton("Data CSV")
        self.rep3_path_button.setMaximumWidth(100)
        self.rep3_path_button.clicked.connect(self.select_rep3_input_file)
        self.rep3_path_button.hide()
        self.rep3_label = QLabel("Replica #3 Data CSV File Path:")
        self.rep3_label.hide()
        self.mut3_label = QLabel("Mutant #2 Data CSV File Path:")
        self.mut3_label.hide()
        layout.addWidget(self.mut3_label, 4, 0)
        layout.addWidget(self.rep3_label, 4, 0)
        layout.addWidget(self.rep3_path_input, 4, 1)
        layout.addWidget(self.rep3_path_button, 4, 2)

        self.rep4_path_input = QLineEdit()
        self.rep4_path_input.setMaximumWidth(225)
        self.rep4_path_input.hide()
        self.rep4_path_button = QPushButton("Data CSV")
        self.rep4_path_button.setMaximumWidth(100)
        self.rep4_path_button.clicked.connect(self.select_rep4_input_file)
        self.rep4_path_button.hide()
        self.rep4_label = QLabel("Replica #4 Data CSV File Path:")
        self.rep4_label.hide()
        self.mut4_label = QLabel("Mutant #3 Data CSV File Path:")
        self.mut4_label.hide()
        layout.addWidget(self.mut4_label, 5, 0)
        layout.addWidget(self.rep4_label, 5, 0)
        layout.addWidget(self.rep4_path_input, 5, 1)
        layout.addWidget(self.rep4_path_button, 5, 2)

        self.rep5_path_input = QLineEdit()
        self.rep5_path_input.setMaximumWidth(225)
        self.rep5_path_input.hide()
        self.rep5_path_button = QPushButton("Data CSV")
        self.rep5_path_button.setMaximumWidth(100)
        self.rep5_path_button.clicked.connect(self.select_rep5_input_file)
        self.rep5_path_button.hide()
        self.rep5_label = QLabel("Replica #5 Data CSV File Path:")
        self.rep5_label.hide()
        self.mut5_label = QLabel("Mutant #4 Data CSV File Path:")
        self.mut5_label.hide()
        layout.addWidget(self.mut5_label, 6, 0)
        layout.addWidget(self.rep5_label, 6, 0)
        layout.addWidget(self.rep5_path_input, 6, 1)
        layout.addWidget(self.rep5_path_button, 6, 2)
       
        self.num_blank_combo = QComboBox()
        self.num_blank_combo.setMaximumWidth(150)
        self.num_blank_combo.addItems(["Select", "1", "2", "3", "4", "5"])
        layout.addWidget(QLabel("Number of Blanks to Subtract:"), 8, 0)
        layout.addWidget(self.num_blank_combo, 8, 1)
        self.num_blank_combo.currentTextChanged.connect(self.toggle_blank_rep_input)
        
        self.blank1_path_input = QLineEdit()
        self.blank1_path_input.setMaximumWidth(225)
        self.blank1_path_input.hide()
        self.blank1_path_button = QPushButton("Blank CSV")
        self.blank1_path_button.setMaximumWidth(100)
        self.blank1_path_button.clicked.connect(self.select_blank1_input_file)
        self.blank1_path_button.hide()
        self.blank1_label = QLabel("Replica #1 Blank CSV File Path:")
        self.blank1_label.hide()
        layout.addWidget(self.blank1_label, 9, 0)
        layout.addWidget(self.blank1_path_input, 9, 1)
        layout.addWidget(self.blank1_path_button, 9, 2)
        
        self.blank2_path_input = QLineEdit()
        self.blank2_path_input.setMaximumWidth(225)
        self.blank2_path_input.hide()
        self.blank2_path_button = QPushButton("Blank CSV")
        self.blank2_path_button.setMaximumWidth(100)
        self.blank2_path_button.clicked.connect(self.select_blank2_input_file)
        self.blank2_path_button.hide()
        self.blank2_label = QLabel("Replica #2 Blank CSV File Path:")
        self.blank2_label.hide()
        layout.addWidget(self.blank2_label, 10, 0)
        layout.addWidget(self.blank2_path_input, 10, 1)
        layout.addWidget(self.blank2_path_button, 10, 2)
        
        self.blank3_path_input = QLineEdit()
        self.blank3_path_input.setMaximumWidth(225)
        self.blank3_path_input.hide()
        self.blank3_path_button = QPushButton("Blank CSV")
        self.blank3_path_button.setMaximumWidth(100)
        self.blank3_path_button.clicked.connect(self.select_blank3_input_file)
        self.blank3_path_button.hide()
        self.blank3_label = QLabel("Replica #3 Blank CSV File Path:")
        self.blank3_label.hide()
        layout.addWidget(self.blank3_label, 11, 0)
        layout.addWidget(self.blank3_path_input, 11, 1)
        layout.addWidget(self.blank3_path_button, 11, 2)
        
        self.blank4_path_input = QLineEdit()
        self.blank4_path_input.setMaximumWidth(225)
        self.blank4_path_input.hide()
        self.blank4_path_button = QPushButton("Blank CSV")
        self.blank4_path_button.setMaximumWidth(100)
        self.blank4_path_button.clicked.connect(self.select_blank4_input_file)
        self.blank4_path_button.hide()
        self.blank4_label = QLabel("Replica #4 Data CSV File Path:")
        self.blank4_label.hide()
        layout.addWidget(self.blank4_label, 12, 0)
        layout.addWidget(self.blank4_path_input, 12, 1)
        layout.addWidget(self.blank4_path_button, 12, 2)

        self.blank5_path_input = QLineEdit()
        self.blank5_path_input.setMaximumWidth(225)
        self.blank5_path_input.hide()
        self.blank5_path_button = QPushButton("Blank CSV")
        self.blank5_path_button.setMaximumWidth(100)
        self.blank5_path_button.clicked.connect(self.select_blank5_input_file)
        self.blank5_path_button.hide()
        self.blank5_label = QLabel("Replica #5 Blank CSV File Path:")
        self.blank5_label.hide()
        layout.addWidget(self.blank5_label, 13, 0)
        layout.addWidget(self.blank5_path_input, 13, 1)
        layout.addWidget(self.blank5_path_button, 13, 2)
        
        
        self.output_dir_input = QLineEdit()
        self.output_dir_input.setMaximumWidth(225)
        self.output_dir_input_button = QPushButton("Output Folder")
        self.output_dir_input_button.setMaximumWidth(100)
        self.output_dir_input_button.clicked.connect(self.select_output_file)
        layout.addWidget(QLabel("Output Directory Path:"), 14, 0)
        layout.addWidget(self.output_dir_input, 14, 1)
        layout.addWidget(self.output_dir_input_button, 14, 2)
        
        self.ok_button = QPushButton("Done Entering Data")
        self.ok_button.clicked.connect(self.ok_button_function)
        layout.addWidget(self.ok_button, 15, 2)
        self.setLayout(layout)

    def toggle_data_rep_input(self, text):
        if text == '2':
            if self.mut_check.isChecked():
                self.wt_label.show()
                self.rep1_label.hide()
                self.rep2_label.hide()
                self.rep3_label.hide()
                self.rep4_label.hide()
                self.rep5_label.hide()
                self.mut2_label.show()
                self.mut3_label.hide()
                self.mut4_label.hide()
                self.mut5_label.hide()
            else:
                self.rep1_label.show()
                self.rep2_label.show()
                self.rep3_label.hide()
                self.rep4_label.hide()
                self.rep5_label.hide()
                self.wt_label.hide()
                self.mut2_label.hide()
                self.mut3_label.hide()
                self.mut4_label.hide()
                self.mut5_label.hide()
            self.rep1_path_input.show()
            self.rep1_path_button.show()
            self.rep2_path_input.show()
            self.rep2_path_button.show()
            self.rep3_path_input.hide()
            self.rep3_path_button.hide()
            self.rep4_path_input.hide()
            self.rep4_path_button.hide()
            self.rep5_path_input.hide()
            self.rep5_path_button.hide()
        elif text == '3':
            if self.mut_check.isChecked():
                self.wt_label.show()
                self.rep1_label.hide()
                self.rep2_label.hide()
                self.rep3_label.hide()
                self.rep4_label.hide()
                self.rep5_label.hide()
                self.mut2_label.show()
                self.mut3_label.show()
                self.mut4_label.hide()
                self.mut5_label.hide()
            else:
                self.rep1_label.show()
                self.rep2_label.show()
                self.rep3_label.show()
                self.rep4_label.hide()
                self.rep5_label.hide()
                self.wt_label.hide()
                self.mut2_label.hide()
                self.mut3_label.hide()
                self.mut4_label.hide()
                self.mut5_label.hide()
            self.rep1_path_input.show()
            self.rep1_path_button.show()
            self.rep2_path_input.show()
            self.rep2_path_button.show()
            self.rep3_path_input.show()
            self.rep3_path_button.show()
            self.rep4_path_input.hide()
            self.rep4_path_button.hide()
            self.rep5_path_input.hide()
            self.rep5_path_button.hide()
        elif text == '4':
            if self.mut_check.isChecked():
                self.wt_label.show()
                self.rep1_label.hide()
                self.rep2_label.hide()
                self.rep3_label.hide()
                self.rep4_label.hide()
                self.rep5_label.hide()
                self.mut2_label.show()
                self.mut3_label.show()
                self.mut4_label.show()
                self.mut5_label.hide()
            else:
                self.rep1_label.show()
                self.rep2_label.show()
                self.rep3_label.show()
                self.rep4_label.show()
                self.rep5_label.hide()
                self.wt_label.hide()
                self.mut2_label.hide()
                self.mut3_label.hide()
                self.mut4_label.hide()
                self.mut5_label.hide()
            self.rep1_path_input.show()
            self.rep1_path_button.show()
            self.rep2_path_input.show()
            self.rep2_path_button.show()
            self.rep3_path_input.show()
            self.rep3_path_button.show()
            self.rep4_path_input.show()
            self.rep4_path_button.show()
            self.rep5_path_input.hide()
            self.rep5_path_button.hide()
        elif text == '5':
            if self.mut_check.isChecked():
                self.wt_label.show()
                self.rep1_label.hide()
                self.rep2_label.hide()
                self.rep3_label.hide()
                self.rep4_label.hide()
                self.rep5_label.hide()
                self.mut2_label.show()
                self.mut3_label.show()
                self.mut4_label.show()
                self.mut5_label.show()
            else:
                self.rep1_label.show()
                self.rep2_label.show()
                self.rep3_label.show()
                self.rep4_label.show()
                self.rep5_label.show()
                self.wt_label.hide()
                self.mut2_label.hide()
                self.mut3_label.hide()
                self.mut4_label.hide()
                self.mut5_label.hide()
            self.rep1_path_input.show()
            self.rep1_path_button.show()
            self.rep2_path_input.show()
            self.rep2_path_button.show()
            self.rep3_path_input.show()
            self.rep3_path_button.show()
            self.rep4_path_input.show()
            self.rep4_path_button.show()
            self.rep5_path_input.show()
            self.rep5_path_button.show()
        else:
            self.rep1_path_input.hide()
            self.rep1_path_button.hide()
            self.rep1_label.hide()
            self.rep2_path_input.hide()
            self.rep2_path_button.hide()
            self.rep2_label.hide()
            self.rep3_path_input.hide()
            self.rep3_path_button.hide()
            self.rep3_label.hide()
            self.rep4_path_input.hide()
            self.rep4_path_button.hide()
            self.rep4_label.hide()
            self.rep5_path_input.hide()
            self.rep5_path_button.hide()
            self.rep5_label.hide()
            self.wt_label.hide()
            self.mut2_label.hide()
            self.mut3_label.hide()
            self.mut4_label.hide()
            self.mut5_label.hide()

    def toggle_blank_rep_input(self, text):
        if text == '1':
            self.blank1_path_input.show()
            self.blank1_path_button.show()
            self.blank1_label.show()
            self.blank2_path_input.hide()
            self.blank2_path_button.hide()
            self.blank2_label.hide()
            self.blank3_path_input.hide()
            self.blank3_path_button.hide()
            self.blank3_label.hide()
            self.blank4_path_input.hide()
            self.blank4_path_button.hide()
            self.blank4_label.hide()
            self.blank5_path_input.hide()
            self.blank5_path_button.hide()
            self.blank5_label.hide()
        elif text == '2':
            self.blank1_path_input.show()
            self.blank1_path_button.show()
            self.blank1_label.show()
            self.blank2_path_input.show()
            self.blank2_path_button.show()
            self.blank2_label.show()
            self.blank3_path_input.hide()
            self.blank3_path_button.hide()
            self.blank3_label.hide()
            self.blank4_path_input.hide()
            self.blank4_path_button.hide()
            self.blank4_label.hide()
            self.blank5_path_input.hide()
            self.blank5_path_button.hide()
            self.blank5_label.hide()
        elif text == '3':
            self.blank1_path_input.show()
            self.blank1_path_button.show()
            self.blank1_label.show()
            self.blank2_path_input.show()
            self.blank2_path_button.show()
            self.blank2_label.show()
            self.blank3_path_input.show()
            self.blank3_path_button.show()
            self.blank3_label.show()
            self.blank4_path_input.hide()
            self.blank4_path_button.hide()
            self.blank4_label.hide()
            self.blank5_path_input.hide()
            self.blank5_path_button.hide()
            self.blank5_label.hide()
        elif text == '4':
            self.blank1_path_input.show()
            self.blank1_path_button.show()
            self.blank1_label.show()
            self.blank2_path_input.show()
            self.blank2_path_button.show()
            self.blank2_label.show()
            self.blank3_path_input.show()
            self.blank3_path_button.show()
            self.blank3_label.show()
            self.blank4_path_input.show()
            self.blank4_path_button.show()
            self.blank4_label.show()
            self.blank5_path_input.hide()
            self.blank5_path_button.hide()
            self.blank5_label.hide()
        elif text == '5':
            self.blank1_path_input.show()
            self.blank1_path_button.show()
            self.blank1_label.show()
            self.blank2_path_input.show()
            self.blank2_path_button.show()
            self.blank2_label.show()
            self.blank3_path_input.show()
            self.blank3_path_button.show()
            self.blank3_label.show()
            self.blank4_path_input.show()
            self.blank4_path_button.show()
            self.blank4_label.show()
            self.blank5_path_input.show()
            self.blank5_path_button.show()
            self.blank5_label.show()
        else:
            self.blank1_path_input.hide()
            self.blank1_path_button.hide()
            self.blank1_label.hide()
            self.blank2_path_input.hide()
            self.blank2_path_button.hide()
            self.blank2_label.hide()
            self.blank3_path_input.hide()
            self.blank3_path_button.hide()
            self.blank3_label.hide()
            self.blank4_path_input.hide()
            self.blank4_path_button.hide()
            self.blank4_label.hide()
            self.blank5_path_input.hide()
            self.blank5_path_button.hide()
            self.blank5_label.hide()

    def select_rep1_input_file(self):
        options = QFileDialog.Options()
        path, _ = QFileDialog.getOpenFileName(self, "Select CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if path is not None:
            self.rep1_path_input.setText(path)
        else:
            QMessageBox.warning(self, 'Critical', 'More Inputs Required')
    def select_rep2_input_file(self):
        options = QFileDialog.Options()
        path, _ = QFileDialog.getOpenFileName(self, "Select CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if path is not None:
            self.rep2_path_input.setText(path)
        else:
            QMessageBox.warning(self, 'Critical', 'More Inputs Required')
    def select_rep3_input_file(self):
        options = QFileDialog.Options()
        path, _ = QFileDialog.getOpenFileName(self, "Select CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if path is not None:
            self.rep3_path_input.setText(path)
        else:
            QMessageBox.warning(self, 'Critical', 'More Inputs Required')
    def select_rep4_input_file(self):
        options = QFileDialog.Options()
        path, _ = QFileDialog.getOpenFileName(self, "Select CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if path is not None:
            self.rep4_path_input.setText(path)
        else:
            QMessageBox.warning(self, 'Critical', 'More Inputs Required')
    def select_rep5_input_file(self):
        options = QFileDialog.Options()
        path, _ = QFileDialog.getOpenFileName(self, "Select CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if path is not None:
            self.rep5_path_input.setText(path)
        else:
            QMessageBox.warning(self, 'Critical', 'More Inputs Required')

    def select_blank1_input_file(self):
        options = QFileDialog.Options()
        path, _ = QFileDialog.getOpenFileName(self, "Select CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if path is not None:
            self.blank1_path_input.setText(path)
        else:
            QMessageBox.warning(self, 'Critical', 'More Inputs Required')
    def select_blank2_input_file(self):
        options = QFileDialog.Options()
        path, _ = QFileDialog.getOpenFileName(self, "Select CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if path is not None:
            self.blank2_path_input.setText(path)
        else:
            QMessageBox.warning(self, 'Critical', 'More Inputs Required')
    def select_blank3_input_file(self):
        options = QFileDialog.Options()
        path, _ = QFileDialog.getOpenFileName(self, "Select CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if path is not None:
            self.blank3_path_input.setText(path)
        else:
            QMessageBox.warning(self, 'Critical', 'More Inputs Required')
    def select_blank4_input_file(self):
        options = QFileDialog.Options()
        path, _ = QFileDialog.getOpenFileName(self, "Select CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if path is not None:
            self.blank4_path_input.setText(path)
        else:
            QMessageBox.warning(self, 'Critical', 'More Inputs Required')
    def select_blank5_input_file(self):
        options = QFileDialog.Options()
        path, _ = QFileDialog.getOpenFileName(self, "Select CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if path is not None:
            self.blank5_path_input.setText(path)
        else:
            QMessageBox.warning(self, 'Critical', 'More Inputs Required')

    def select_output_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        directory = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if directory:
            self.output_dir_input.setText(directory)
            self.folder_selected.emit(directory)
            os.environ['WORKING_DIR'] = directory
        else:
            QMessageBox.warning(self, 'Critical', 'More Inputs Required')

    def mut_toggle(self):
        if self.mut_check.isChecked():
            text = self.num_rep_combo.currentText()
            if text == '2':
                self.wt_label.show()
                self.rep1_label.hide()
                self.rep2_label.hide()
                self.rep3_label.hide()
                self.rep4_label.hide()
                self.rep5_label.hide()
                self.mut2_label.show()
                self.mut3_label.hide()
                self.mut4_label.hide()
                self.mut5_label.hide()
            elif text == '3':
                self.wt_label.show()
                self.rep1_label.hide()
                self.rep2_label.hide()
                self.rep3_label.hide()
                self.rep4_label.hide()
                self.rep5_label.hide()
                self.mut2_label.show()
                self.mut3_label.show()
                self.mut4_label.hide()
                self.mut5_label.hide()
            elif text == '4':
                self.wt_label.show()
                self.rep1_label.hide()
                self.rep2_label.hide()
                self.rep3_label.hide()
                self.rep4_label.hide()
                self.rep5_label.hide()
                self.mut2_label.show()
                self.mut3_label.show()
                self.mut4_label.show()
                self.mut5_label.hide()
            elif text == '5':
                self.wt_label.show()
                self.rep1_label.hide()
                self.rep2_label.hide()
                self.rep3_label.hide()
                self.rep4_label.hide()
                self.rep5_label.hide()
                self.mut2_label.show()
                self.mut3_label.show()
                self.mut4_label.show()
                self.mut5_label.show()
            else:
                self.wt_label.hide()
                self.rep1_label.hide()
                self.rep2_label.hide()
                self.rep3_label.hide()
                self.rep4_label.hide()
                self.rep5_label.hide()
                self.mut2_label.hide()
                self.mut3_label.hide()
                self.mut4_label.hide()
                self.mut5_label.hide()
        else:
            text = self.num_rep_combo.currentText()
            if text == '2':
                self.wt_label.hide()
                self.rep1_label.show()
                self.rep2_label.show()
                self.rep3_label.hide()
                self.rep4_label.hide()
                self.rep5_label.hide()
                self.mut2_label.hide()
                self.mut3_label.hide()
                self.mut4_label.hide()
                self.mut5_label.hide()
            elif text == '3':
                self.wt_label.hide()
                self.rep1_label.show()
                self.rep2_label.show()
                self.rep3_label.show()
                self.rep4_label.hide()
                self.rep5_label.hide()
                self.mut2_label.hide()
                self.mut3_label.hide()
                self.mut4_label.hide()
                self.mut5_label.hide()
            elif text == '4':
                self.wt_label.hide()
                self.rep1_label.show()
                self.rep2_label.show()
                self.rep3_label.show()
                self.rep4_label.show()
                self.rep5_label.hide()
                self.mut2_label.hide()
                self.mut3_label.hide()
                self.mut4_label.hide()
                self.mut5_label.hide()
            elif text == '5':
                self.wt_label.hide()
                self.rep1_label.show()
                self.rep2_label.show()
                self.rep3_label.show()
                self.rep4_label.show()
                self.rep5_label.show()
                self.mut2_label.hide()
                self.mut3_label.hide()
                self.mut4_label.hide()
                self.mut5_label.hide()
            else:
                self.wt_label.hide()
                self.rep1_label.hide()
                self.rep2_label.hide()
                self.rep3_label.hide()
                self.rep4_label.hide()
                self.rep5_label.hide()
                self.mut2_label.hide()
                self.mut3_label.hide()
                self.mut4_label.hide()
                self.mut5_label.hide()
 

    def ok_button_function(self):
        if self.rep1_path_input.text().strip():
            pathd = os.path.join(self.output_dir_input.text(), "path_rep_data.txt")
            with open(pathd, "w") as f:
                f.write(self.rep1_path_input.text())
                f.write("\n")
        if self.rep2_path_input.text().strip():
            pathd = os.path.join(self.output_dir_input.text(), "path_rep_data.txt")
            with open(pathd, "a") as f:
                f.write(self.rep2_path_input.text())
                f.write("\n")
        if self.rep3_path_input.text().strip():
            pathd = os.path.join(self.output_dir_input.text(), "path_rep_data.txt")
            with open(pathd, "a") as f:
                f.write(self.rep3_path_input.text())
                f.write("\n")
        if self.rep4_path_input.text().strip():
            pathd = os.path.join(self.output_dir_input.text(), "path_rep_data.txt")
            with open(pathd, "a") as f:
                f.write(self.rep4_path_input.text())
                f.write("\n")
        if self.rep5_path_input.text().strip():
            pathd = os.path.join(self.output_dir_input.text(), "path_rep_data.txt")
            with open(pathd, "a") as f:
                f.write(self.rep5_path_input.text())
                f.write("\n")
        if self.blank1_path_input.text().strip():
            pathb = os.path.join(self.output_dir_input.text(), "blank_rep_data.txt")
            with open(pathb, "w") as f:
                f.write(self.blank1_path_input.text())
                f.write("\n")
        if self.blank2_path_input.text().strip():
            pathd = os.path.join(self.output_dir_input.text(), "blank_rep_data.txt")
            with open(pathb, "a") as f:
                f.write(self.blank2_path_input.text())
                f.write("\n")
        if self.blank3_path_input.text().strip():
            pathd = os.path.join(self.output_dir_input.text(), "blank_rep_data.txt")
            with open(pathb, "a") as f:
                f.write(self.blank3_path_input.text())
                f.write("\n")
        if self.blank4_path_input.text().strip():
            pathd = os.path.join(self.output_dir_input.text(), "blank_rep_data.txt")
            with open(pathb, "a") as f:
                f.write(self.blank4_path_input.text())
                f.write("\n")
        if self.blank5_path_input.text().strip():
            pathd = os.path.join(self.output_dir_input.text(), "blank_rep_data.txt")
            with open(pathb, "a") as f:
                f.write(self.blank5_path_input.text())
                f.write("\n")
        
        if self.mut_check.isChecked():
            mut = os.path.join(self.output_dir_input.text(), "mutant.txt")
            with open(mut, "w") as f:
                f.write("mutant\n")
        else:
            mut = os.path.join(self.output_dir_input.text(), "mutant.txt")
            with open(mut, "w") as f:
                f.write("\n") 

        self.done.emit()
        self.close()

class KineticAnalysisTool(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Kinetic Analysis Tool (KAT)")
        self.setGeometry(100, 100, 1200, 900)

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        palette = self.central_widget.palette()
        palette.setColor(QPalette.Window, QColor("#f0f4f8"))
        self.central_widget.setAutoFillBackground(True)
        self.central_widget.setPalette(palette)

        # Main horizontal layout: left (inputs), middle (columns + instructions), right (tabs)
        main_layout = QHBoxLayout(self.central_widget)

        # Left side layout: inputs
        left_layout = QVBoxLayout()

        input_form = QGridLayout()

        self.output_name_input = QLineEdit()
        self.output_name_input.setMaximumWidth(150)
        input_form.addWidget(QLabel("Output Filename:"), 0, 0)
        input_form.addWidget(self.output_name_input, 0, 1)
        self.batch_input = QPushButton("Replica Inputs")
        self.batch_input.setStyleSheet("background-color: #f4c178; color: black;")
        self.batch_input.setMaximumWidth(100)
        self.batch_input.clicked.connect(self.batch_input_window)
        input_form.addWidget(self.batch_input, 0, 2)

        self.file_path_input = QLineEdit()
        self.file_path_input.setMaximumWidth(225)
        self.file_path_button = QPushButton("Data CSV")
        self.file_path_button.setMaximumWidth(100)
        self.file_path_button.clicked.connect(self.select_input_file)
        input_form.addWidget(QLabel("Data CSV File Path:"), 1, 0)
        input_form.addWidget(self.file_path_input, 1, 1)
        input_form.addWidget(self.file_path_button, 1, 2)

        self.blank_path_input = QLineEdit()
        self.blank_path_input.setMaximumWidth(225)
        self.blank_path_input.setPlaceholderText("only if applicable")
        self.blank_path_button = QPushButton("Blank CSV")
        self.blank_path_button.setMaximumWidth(100)
        self.blank_path_button.clicked.connect(self.select_blank_file)
        input_form.addWidget(QLabel("Blank CSV File Path:"), 2, 0)
        input_form.addWidget(self.blank_path_input, 2, 1)
        input_form.addWidget(self.blank_path_button, 2, 2)

        self.output_dir_input = QLineEdit()
        self.output_dir_input.setMaximumWidth(225)
        self.output_dir_input_button = QPushButton("Output Folder")
        self.output_dir_input_button.setMaximumWidth(100)
        self.output_dir_input_button.clicked.connect(self.select_output_file)
        input_form.addWidget(QLabel("Output Directory Path:"), 3, 0)
        input_form.addWidget(self.output_dir_input, 3, 1)
        input_form.addWidget(self.output_dir_input_button, 3, 2)

        self.data_type_combo = QComboBox()
        self.data_type_combo.setMaximumWidth(150)
        self.data_type_combo.addItems(["Fluorescence", "Absorbance"])
        input_form.addWidget(QLabel("Kinetic Data Type:"), 4, 0)
        input_form.addWidget(self.data_type_combo, 4, 1)

        self.absorbance_input = QLineEdit()
        self.absorbance_input.setPlaceholderText("Enter absorbance details")
        self.absorbance_input.setValidator(QDoubleValidator())
        self.absorbance_input.hide()  # Initially hidden
        self.absorbance_label = QLabel("Molar Absorptivity:")
        self.absorbance_label.hide()
        input_form.addWidget(self.absorbance_label, 5, 0)
        input_form.addWidget(self.absorbance_input, 5, 1)

        # Connect signal to handler
        self.data_type_combo.currentTextChanged.connect(self.toggle_absorbance_input)

        self.enzyme_conc_input = QLineEdit()
        self.enzyme_conc_input.setMaximumWidth(150)
        self.enzyme_conc_input.setValidator(QDoubleValidator())
        input_form.addWidget(QLabel("Enzyme Concentration:"), 6, 0)
        input_form.addWidget(self.enzyme_conc_input, 6, 1)

        left_layout.addLayout(input_form)

        # Substrate Info Group
        substrate_group = QGroupBox("Substrate Information")
        substrate_group.setStyleSheet("QGroupBox { background-color: #b4dfff; border: 1px solid #4db2ff; padding: 10px; }")
        substrate_group.setMaximumHeight(120)
        substrate_layout = QGridLayout()
        self.substrate_number = QLineEdit()
        self.substrate_number.setMaximumWidth(100)
        self.substrate_number.setValidator(QDoubleValidator())
        self.substrate_dilution = QLineEdit()
        self.substrate_dilution.setMaximumWidth(100)
        self.substrate_dilution.setValidator(QDoubleValidator())
        self.substrate_max = QLineEdit()
        self.substrate_max.setMaximumWidth(100)
        self.substrate_max.setValidator(QDoubleValidator())
        
        substrate_layout.addWidget(QLabel("# Substrate Concentrations:"), 0, 0)
        substrate_layout.addWidget(self.substrate_number, 0, 1)
        substrate_layout.addWidget(QLabel("Dilution Factor:"), 1, 0)
        substrate_layout.addWidget(self.substrate_dilution, 1, 1)
        substrate_layout.addWidget(QLabel("Max Substrate Concentration:"), 2, 0)
        substrate_layout.addWidget(self.substrate_max, 2, 1)
        substrate_group.setLayout(substrate_layout)
        left_layout.addWidget(substrate_group)

        # Time-course Info Group
        time_group = QGroupBox("Time-Course Information")
        time_group.setStyleSheet("QGroupBox { background-color: #ffc5c0; border: 1px solid #ff746a; padding: 10px; }")
        time_layout = QGridLayout()
        
        time_help_icon = QApplication.style().standardIcon(QStyle.SP_MessageBoxInformation)
        self.time_help_button = QPushButton()
        self.time_help_button.setFixedSize(20, 15)
        self.time_help_button.setIcon(time_help_icon)
        self.time_help_button.setIconSize(QSize(24, 24))
        self.time_help_button.setToolTip("Help with Inputting Time Information")
        self.time_help_button.clicked.connect(self.open_time_help_window)
        self.time_help_layout = QHBoxLayout()
        self.time_help_layout.addStretch()
        self.time_help_layout.addWidget(self.time_help_button)
        time_layout.addLayout(self.time_help_layout, 0, 3)
        
        self.auto_lin_range = QCheckBox("Auto Calculate Linear Range")
        self.auto_lin_range.setChecked(True)  # Optional: set initial state
        self.auto_lin_range.stateChanged.connect(self.on_checkbox_toggled)        
        
        self.v_win = QLineEdit()
        self.v_win.setMaximumWidth(100)
        self.v_win.setValidator(QDoubleValidator())
        self.time_min = QLineEdit()
        self.time_min.setMaximumWidth(100)
        self.time_min.setValidator(QDoubleValidator())
        self.time_min.hide()
        self.time_max = QLineEdit()
        self.time_max.setMaximumWidth(100)
        self.time_max.setValidator(QDoubleValidator())
        self.time_max.hide()
        self.time_step = QLineEdit()
        self.time_step.setMaximumWidth(100)
        self.time_step.setValidator(QDoubleValidator())
        self.time_step.hide()

        time_layout.addWidget(self.auto_lin_range, 0, 0)
        self.v_win_label = QLabel("Computed Velocity Window:")
        time_layout.addWidget(self.v_win_label, 1, 0)
        time_layout.addWidget(self.v_win, 1, 1)
        self.time_min_label = QLabel("Time Min (min):")
        self.time_min_label.hide()
        time_layout.addWidget(self.time_min_label, 2, 0)
        time_layout.addWidget(self.time_min, 2, 1)
        self.time_max_label = QLabel("Time Max (min:")
        self.time_max_label.hide()
        time_layout.addWidget(self.time_max_label, 3, 0)
        time_layout.addWidget(self.time_max, 3, 1)
        self.time_step_label = QLabel("Time Step (min):")
        self.time_step_label.hide()
        time_layout.addWidget(self.time_step_label, 4, 0)
        time_layout.addWidget(self.time_step, 4, 1)
        time_group.setLayout(time_layout)
        left_layout.addWidget(time_group)

        logos_group = QGroupBox()
        logos_group.setStyleSheet("QGroupBox { background-color: #f0f4f8; border: 1px solid #f0f4f8; padding: 10px; }")
        logos_layout = QGridLayout()
        pixmap_logo = QPixmap(get_resource_path(".imgs/logo_v3.png"))
        self.logo_label = QLabel()
        self.logo_label.setPixmap(pixmap_logo)
        self.logo_label.setScaledContents(True)
        self.logo_label.setFixedSize(143, 130)
        self.opacity_effect = QGraphicsOpacityEffect(self.logo_label)
        self.opacity_effect.setOpacity(0.4)  # Set opacity (0.0 for fully transparent, 1.0 for fully opaque)
        self.logo_label.setGraphicsEffect(self.opacity_effect)
        logos_layout.addWidget(self.logo_label, 0, 2)
        logos_group.setLayout(logos_layout)
        logos_group.setMaximumHeight(175)
        left_layout.addWidget(logos_group)

        main_layout.addLayout(left_layout) 

        # Middle layout: Column Selector + Instructions
        middle_layout = QVBoxLayout()

        ascii_art = """
    ╭────────────────────────────────────────╮
    │  ██╗  ██╗ █████╗ ████████╗             │
    │  ██║ ██╔╝██╔══██╗╚══██╔══╝             │
    │  █████╔╝ ███████║   ██║    █████╗      │
    │  ██╔═██╗ ██╔══██║   ██║    ╚════╝      │
    │  ██║  ██╗██║  ██║   ██║                │
    │  ╚═╝  ╚═╝╚═╝  ╚═╝   ╚═╝                │
    │      KINETIC ANALYSIS TOOLKIT          │
    ╰────────────────────────────────────────╯
        """

        title_logo = QLabel(ascii_art)
        title_logo.setAlignment(Qt.AlignCenter)
        title_logo.setStyleSheet("""
            QLabel {
            color: #23395d;
            background-color: "#f0f4f8";
            font-family: 'Courier New', monospace;
            font-size: 10pt;
            padding: 8px;
            border: 3px solid #23395d;
            border-radius: 10px;
            }
            """)
        title_logo.setMaximumHeight(150)
        title_logo.setMinimumWidth(670)
        
        pixmap_title = QPixmap(get_resource_path(".imgs/title_logo.png")).scaled(470, 125, Qt.IgnoreAspectRatio, Qt.SmoothTransformation)
        # Create a QLabel and set the pixmap
        title_label = QLabel()
        title_label.setPixmap(pixmap_title)
        title_label.setAlignment(Qt.AlignCenter)

        # Optionally set fixed size for the label
        title_label.setFixedSize(670, 150)  # adjust width, height to your liking

        # Add this label to a layout, e.g. next to your tabs or buttons

        middle_layout.addWidget(title_logo)

        self.output_name_input.setFixedHeight(25)
        self.file_path_input.setFixedHeight(25)
        self.substrate_number.setFixedHeight(22)
        self.substrate_dilution.setFixedHeight(22)
        self.substrate_max.setFixedHeight(22)
        self.time_min.setFixedHeight(20)
        self.time_max.setFixedHeight(20)
        self.time_step.setFixedHeight(20)
        self.data_type_combo.setFixedHeight(25)
        self.absorbance_input.setFixedHeight(25)
        self.enzyme_conc_input.setFixedHeight(25)

        # Show output values
        output_data_group = QGroupBox("Kinetic Parameters")
        output_data_group.setStyleSheet("QGroupBox { background-color: #f0f4f8; border: 1px solid #94cc4c; padding: 10px; }")
        output_data_layout = QGridLayout()
        
        self.h_output = QLineEdit()
        self.h_output.setMaximumWidth(150)
        output_data_layout.addWidget(QLabel("Hill Coefficient:"), 0, 0)
        output_data_layout.addWidget(self.h_output, 1, 0)

        self.vmax_output = QLineEdit()
        self.output_name_input.setMaximumWidth(150)
        output_data_layout.addWidget(QLabel("Vmax Value (/min):"), 0, 1)
        output_data_layout.addWidget(self.vmax_output, 1, 1)

        self.km_output = QLineEdit()
        self.km_output.setMaximumWidth(150)
        output_data_layout.addWidget(QLabel("Km Value:"), 0, 2)
        output_data_layout.addWidget(self.km_output, 1, 2)

        self.kcat_output = QLineEdit()
        self.kcat_output.setMaximumWidth(150)
        output_data_layout.addWidget(QLabel("Kcat Value (/sec):"), 0, 3)
        output_data_layout.addWidget(self.kcat_output, 1, 3)

        output_data_group.setLayout(output_data_layout)
        output_data_group.setMaximumHeight(90)
        middle_layout.addWidget(output_data_group)

        # Add graph display button below instructions (if you want it here)
        self.graph_button = QPushButton("Display Graph and Output Values")
        self.graph_button.setStyleSheet("background-color: #3399cc; color: white;")
        middle_layout.addWidget(self.graph_button)
        self.graph_button.clicked.connect(self.display_graph)

        # Image Display Area (Center)
        
        self.display_group = QGroupBox("Graph Display")
        self.display_group.setStyleSheet("QGroupBox { background-color: #ffffff; border: 1px solid #cccccc; padding: 10px; }")
        self.display_layout = QGridLayout()
        
        self.image_display = QLabel()
        self.image_display.setStyleSheet("background-color: white; border: 1px solid #ccc;")
        self.image_display.setMinimumHeight(400)
        self.image_display.setMinimumWidth(670)

        middle_layout.insertWidget(1, self.image_display)

        main_layout.addLayout(middle_layout)

        # Right side layout: Tabs
        right_layout = QVBoxLayout()

        self.tabs = QTabWidget()
        self.tabs.setStyleSheet("QTabWidget::pane { background: #dfd4ff; border: 1px solid #b498ff; padding: 10px }")
        self.tabs.setMaximumHeight(400)
        self.tab_classical = QWidget()
        self.tab_complex = QWidget()
        self.tab_inhib = QWidget()

        self.tabs.addTab(self.tab_classical, "Classical Models")
        self.tabs.addTab(self.tab_complex, "Complex Models")
        """self.tabs.addTab(self.tab_inhib, "Inhibition")"""

        self.setup_tabs()
        right_layout.addWidget(self.tabs)


        pixmap = QPixmap(get_resource_path(".imgs/Cat Scientist at Work.png"))

        # Create a QLabel and set the pixmap
        image_label = QLabel()
        image_label.setPixmap(pixmap)
        image_label.setAlignment(Qt.AlignCenter)
        image_label.setScaledContents(True)  # optional: scale image to fit label size

        # Optionally set fixed size for the label
        image_label.setFixedSize(325, 450)  # adjust width, height to your liking

        # Add this label to a layout, e.g. next to your tabs or buttons
        right_layout.addWidget(image_label)
        
        main_layout.addLayout(right_layout)

        self.help_button = QPushButton("Help")
        self.help_button.setStyleSheet("background-color: #C41E3A; color: white;")
        self.help_button.clicked.connect(self.open_gen_help_window)
        right_layout.addWidget(self.help_button)
        

        # Welcome label centered at bottom spanning width
        self.welcome_label = QLabel(
            "Welcome to KAT, the Kinetic Analysis Toolkit!\n"
            "Load CSV, set parameters, then run analysis."
        )
        self.welcome_label.setFont(QFont("Times New Roman", 20))
        self.welcome_label.setAlignment(Qt.AlignCenter)
        self.welcome_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)

        # Add welcome label to a new vertical layout below main layout
        outer_layout = QVBoxLayout()
        outer_layout.addLayout(main_layout)
        outer_layout.addWidget(self.welcome_label)

        self.central_widget.setLayout(outer_layout)

        self.statusBar().showMessage("Ready")
	
	# Create a label to hold the animation
        self.loading_label = QLabel(self)
        self.loading_label.setFixedSize(175, 175)  # adjust size as needed
        #self.loading_label.setStyleSheet("background: #ffffff;")  # optional
        # Set up the movie (GIF)
        gif_path = get_movie_path('.imgs/kat_loading.gif')
        self.movie = QMovie(gif_path)
        if not self.movie.isValid():
            print(f"QMovie failed to load: {gif_path}")
        else:
            print("QMovie successfully loaded")
        self.movie.setScaledSize(QSize(175, 175))
        self.loading_label.setMovie(self.movie)
        self.loading_label.hide()  # Initially hidden
        self.center_loading_label()

    def center_loading_label(self):
        window_size = self.size()
        label_size = self.loading_label.size()
        x = (window_size.width() - label_size.width()) // 2
        y = (window_size.height() - label_size.height()) // 2
        self.loading_label.move(x, y)
    
    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.center_loading_label()

    def setup_tabs(self):
        classical_layout = QVBoxLayout()
        class_tabs = QTabWidget()
        class_tabs.setStyleSheet("QTabWidget::pane { background: #dfd4ff; border: 1px solid #b498ff; padding: 10px }")
        tab_class_standard = QWidget()
        tab_class_rep = QWidget()
        class_tabs.addTab(tab_class_standard, "Standard")
        class_tabs.addTab(tab_class_rep, "Replicas")

        class_rep_layout = QVBoxLayout()
        class_std_layout = QVBoxLayout()
        hill_averages = os.path.join(SCRIPTS_DIR, 'hill-averages_winset.py')
        class_std_layout.addWidget(self.make_script_button("Hill Kinetics", hill_averages))
        class_std_layout.addWidget(self.make_script_button("Michaelis-Menten Kinetics", os.path.join(SCRIPTS_DIR, 'mm-averages.py')))
        class_rep_layout.addWidget(self.make_script_button("Hill Replicas", os.path.join(SCRIPTS_DIR, "hill-replicas.py")))
        class_rep_layout.addWidget(self.make_script_button("Michaelis-Menten Replicas", os.path.join(SCRIPTS_DIR, "mm-replicas.py")))
        tab_class_standard.setLayout(class_std_layout)
        tab_class_rep.setLayout(class_rep_layout)
        class_tabs.setMaximumWidth(300)
        classical_layout.addWidget(class_tabs)
        self.tab_classical.setLayout(classical_layout)

        complex_layout = QVBoxLayout()
        comp_tabs = QTabWidget()
        comp_tabs.setStyleSheet("QTabWidget::pane { background: #dfd4ff; border: 1px solid #b498ff; padding: 10px }")
        tab_comp_standard = QWidget()
        tab_comp_rep = QWidget()
        comp_tabs.addTab(tab_comp_standard, "Standard")
        comp_tabs.addTab(tab_comp_rep, "Replicas")

        comp_std_layout = QVBoxLayout()
        comp_rep_layout = QVBoxLayout()
        comp_std_layout.addWidget(self.make_script_button("Monod-Wyman-Changeux", os.path.join(SCRIPTS_DIR, 'mwc-averages.py')))
        knf_run = QPushButton("Koshland-Nemethy-Filmer")
        comp_std_layout.addWidget(knf_run)
        self.site_data = QLineEdit()
        self.site_data.setPlaceholderText("enter # of active sites")
        self.site_data.hide()
        self.run_knf_btn = QPushButton("Run KNF")
        self.run_knf_btn.clicked.connect(lambda: self.run_script(os.path.join(SCRIPTS_DIR, 'knf-averages.py')))
        self.run_knf_btn.hide()
        knf_run.clicked.connect(self.show_site_data)
        comp_std_layout.addWidget(self.site_data)
        comp_std_layout.addWidget(self.run_knf_btn)
        tab_comp_standard.setLayout(comp_std_layout)

        comp_rep_layout.addWidget(self.make_script_button("Monod-Wyman-Changeux Replicas", os.path.join(SCRIPTS_DIR, "mwc-replicas.py")))
        knf_run_rep = QPushButton("Koshland-Nemethy-Filmer Replicas")
        comp_rep_layout.addWidget(knf_run_rep)
        self.site_data_rep = QLineEdit()
        self.site_data_rep.setPlaceholderText("enter # of active sites, comma deliminated if necessary")
        self.site_data_rep.hide()
        self.run_knf_btn_rep = QPushButton("Run KNF Replicas")
        self.run_knf_btn_rep.clicked.connect(lambda: self.run_script(os.path.join(SCRIPTS_DIR, 'knf-replicas.py')))
        self.run_knf_btn_rep.hide()
        knf_run_rep.clicked.connect(self.show_site_data_rep)
        comp_rep_layout.addWidget(self.site_data_rep)
        comp_rep_layout.addWidget(self.run_knf_btn_rep)
        tab_comp_rep.setLayout(comp_rep_layout)
        complex_layout.addWidget(comp_tabs)
        
        self.tab_complex.setLayout(complex_layout)
    """
        inhib_layout = QVBoxLayout()
        inhib_layout.addWidget(self.make_script_button("Competitive?", "../scripts/comp-averages.py"))
        inhib_layout.addWidget(self.make_script_button("Uncompetitive?", "../scripts/hill-averages_winset.py"))
        inhib_layout.addWidget(self.make_script_button("Non-Competitive?", "../scripts/hill-averages_winset.py"))
        self.tab_inhib.setLayout(inhib_layout)
    """

    def show_site_data(self):
        self.site_data.show()
        self.run_knf_btn.show()
    
    def show_site_data_rep(self):
        self.site_data_rep.show()
        self.run_knf_btn_rep.show()

    def make_script_button(self, label, script):
        btn = QPushButton(label)
        btn.clicked.connect(lambda: self.run_script(get_resource_path(script)))
        return btn

    def toggle_column(self, button, index):
        if index in self.selected_columns:
            self.selected_columns.remove(index)
            button.setStyleSheet("")
        else:
            self.selected_columns.append(index)
            button.setStyleSheet("background-color: lightgreen")
        self.columns_label.setText(f"Selected Columns: {self.selected_columns}")

    def select_input_file(self):
        options = QFileDialog.Options()
        path, _ = QFileDialog.getOpenFileName(self, "Select CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if path is not None:
            self.file_path_input.setText(path)
        else:
            QMessageBox.warning(self, 'Critical', 'More Inputs Required')
    
    def select_blank_file(self):
        options = QFileDialog.Options()
        path, _ = QFileDialog.getOpenFileName(self, "Select CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if path is not None:
            self.blank_path_input.setText(path)
        else:
            QMessageBox.warning(self, 'Critical', 'More Inputs Required')

    def select_output_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog  
        directory = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if directory:
            self.output_dir_input.setText(directory)
            os.environ['WORKING_DIR'] = directory
        else:
            QMessageBox.warning(self, 'Critical', 'More Inputs Required')
    
    def run_script(self, script_path):
        self.statusBar().showMessage("Running analysis script...")
        self.loading_label.show()
        self.movie.start()
        if self.file_path_input.text() == 'Batch Import':
            QMessageBox.warning(self, 'Warning', 'You have entered replica data, make sure you select a model under the Replicas tab')

        if self.substrate_number.text().strip():
            substrate = os.path.join(self.output_dir_input.text(), "substrate_data.txt")
            with open(substrate, "w") as f:
                f.write(self.substrate_number.text())
                f.write("\n")
        else:
            QMessageBox.warning(self, 'Warning', 'Substrate Number Required')
        if self.substrate_dilution.text().strip():
            substrate = os.path.join(self.output_dir_input.text(), "substrate_data.txt")
            with open(substrate, "a") as f:
                f.write(self.substrate_dilution.text())
                f.write("\n")
        else:
            QMessageBox.warning(self, 'Warning', 'Substrate Dilution Factor Required')
        if self.substrate_max.text().strip():
            substrate = os.path.join(self.output_dir_input.text(), "substrate_data.txt")
            with open(substrate, "a") as f:
                f.write(self.substrate_max.text())
                f.write("\n")
        else:
            QMessageBox.warning(self, 'Warning', 'Substrate Maxiumum Required')
        
        if self.auto_lin_range.isChecked():
            time = os.path.join(self.output_dir_input.text(), "time_data.txt")
            with open(time, "w") as f:
                f.write('True\n')
                f.write(self.v_win.text())
                f.write("\n")
        else:
            time = os.path.join(self.output_dir_input.text(), "time_data.txt")
            with open(time, "w") as f:
                f.write('False\n')
                f.write(self.v_win.text())
                f.write("\n")
                if self.time_min.text().strip():
                    f.write(self.time_min.text())
                    f.write("\n")
                else:
                    QMessageBox.warning(self, 'Warning', 'Time Minimum Required')
                if self.time_max.text().strip():
                    f.write(self.time_max.text())
                    f.write("\n")
                else:
                    QMessageBox.warning(self, 'Warning', 'Time Maximum Required')
                if self.time_step.text().strip():
                    f.write(self.time_step.text())
                    f.write("\n")
                else:
                    QMessageBox.warning(self, 'Warning', 'Time Step Required')

        if self.output_name_input.text().strip():
            name = os.path.join(self.output_dir_input.text(), "name_data.txt")
            with open(name, "w") as f:
                f.write(self.output_name_input.text())
                f.write("\n")
        else:
            QMessageBox.warning(self, 'Warning', 'Name of Output Required')

        if self.data_type_combo.currentText() == "Absorbance":
            data = os.path.join(self.output_dir_input.text(), "data_type.txt")
            with open(data, "w") as f:
                f.write(self.data_type_combo.currentText())
                f.write("\n")
        elif self.data_type_combo.currentText() == "Fluorescence":
            data = os.path.join(self.output_dir_input.text(), "data_type.txt")
            with open(data, "w") as f:
                f.write(self.data_type_combo.currentText())
                f.write("\n")
        else:
            QMessageBox.warning(self, 'Warning', 'Data-Type Required')
        if self.absorbance_input.text().strip():
            data = os.path.join(self.output_dir_input.text(), "data_type.txt")
            with open(data, "a") as f:
                f.write(self.absorbance_input.text())
                f.write("\n")
        if self.file_path_input.text().strip():
            pathd = os.path.join(self.output_dir_input.text(), "path_data.txt")
            with open(pathd, "w") as f:
                f.write(self.file_path_input.text())
                f.write("\n")
                f.write(self.blank_path_input.text())
                f.write("\n")
        else:
            QMessageBox.warning(self, 'Warning', 'Data Path Required')
       
        if self.site_data.text().strip():
            site = os.path.join(self.output_dir_input.text(), "sites.txt")
            with open(site, "w") as f:
                f.write(self.site_data.text())
                f.write("\n")
        if self.site_data_rep.text().strip():
            site = os.path.join(self.output_dir_input.text(), "sites.txt")
            with open(site, "w") as f:
                f.write(self.site_data_rep.text())
                f.write("\n")
        self.site_data.hide()
        self.run_knf_btn.hide()
        self.site_data_rep.hide()
        self.run_knf_btn_rep.hide()

        current_index = self.tabs.currentIndex()
        current_tab_text = self.tabs.tabText(current_index)
        if current_tab_text == "Complex Models" and int(self.substrate_number.text().strip()) < 30:
            QMessageBox.warning(self, "Warning", 'Warning: >30 substrate values is recommended to avoid overfitting complex models')


        self.thread = QThread()
        self.worker = ScriptWorker(script_path)
        self.worker.moveToThread(self.thread)

        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)

        self.worker.finished.connect(self.movie.stop)
        self.worker.finished.connect(self.loading_label.hide)
        self.worker.finished.connect(lambda: self.statusBar().showMessage("Done"))

        self.thread.start()
        
    def do_work(self, script_path):
        try:
            run_internal_script(script_path, call_main=True)
        finally:
            self.movie.stop()
            self.loading_label.hide()

    def display_graph(self):
        filename = self.output_name_input.text().strip()
        output_dir = self.output_dir_input.text().strip()
        if not filename:
            print("No filename specified.")
            return
        if not output_dir:
            output_dir = os.getcwd()

        image_path = os.path.join(output_dir, f"{filename}.png")
        if not os.path.exists(image_path):
            QMessageBox.warning(self, "Error Finding Path", f"{image_path} not found.")
            return

        pixmap = QPixmap(image_path)
        scaled_pixmap = pixmap.scaled(
            self.image_display.width(), self.image_display.height(),
            Qt.KeepAspectRatio, Qt.SmoothTransformation
        )

        self.image_display.setPixmap(scaled_pixmap)

        kp_path = os.path.join(output_dir, f'{filename}.txt')
        if not os.path.exists(kp_path):
            QMessageBox.warning(self, "Error Finding Path", f"{kp_path} not found.")
            return

        if os.path.exists(os.path.join(output_dir, "mutant.txt")):
            with open(os.path.join(output_dir, "mutant.txt"), 'r') as file:
                mut = [line.strip() for line in file.readlines()]
        else: 
            mut = []

        with open(kp_path, "r") as file:
                lines = [line.strip() for line in file.readlines()]
                line_1 = str(lines[0])
                vmax = str(lines[1])
                line_2 = lines[1].split()
                line_3 = str(lines[2])
                output_data = [line_1, vmax, line_3]
        if 'mutant' in mut:
            self.dat = DisplayMutData()
            self.dat.show()
        elif 'Complex model used' in lines:
            self.dat = DisplayComplexData()
            self.dat.show()
        else:
            self.h_output.setText(str(output_data[0]))
            self.vmax_output.setText(str(output_data[1]))
            self.km_output.setText(str(output_data[2]))
            enzyme_conc = self.enzyme_conc_input.text().strip()
            if enzyme_conc:
                kcat_avg = float(line_2[0]) / float(enzyme_conc) / 60
                if len(line_2) > 1:
                    kcat_std = float(line_2[2]) / float(enzyme_conc) / 60
                    val = '%.3f'%(kcat_avg) + '\u00B1' + '%.3f'%(kcat_std)
                    self.kcat_output.setText(str(val))
                else:
                    kcat = float(output_data[1]) / float(enzyme_conc) / 60
                    self.kcat_output.setText(str(kcat))
            else:
                self.kcat_output.setText("N/A")

    def open_gen_help_window(self):
        self.help_window = GeneralHelpWindow()
        self.help_window.show()
    def open_time_help_window(self):
        self.help_window = TimeHelpWindow()
        self.help_window.show()
    def batch_input_window(self):
        self.batch_import = BatchImportCSV()
        self.batch_import.folder_selected.connect(self.set_main_output_folder)
        self.batch_import.done.connect(self.set_batch_import_label)
        self.batch_import.show()
    
    def set_main_output_folder(self, path):
        self.output_dir_input.setText(path)
        
    def set_batch_import_label(self):
        self.file_path_input.setText('Batch Import')
        self.blank_path_input.setText('Batch Import')

    def toggle_absorbance_input(self, text):
        if text == "Absorbance":
            self.absorbance_label.show()
            self.absorbance_input.show()
        else:
            self.absorbance_label.hide()
            self.absorbance_input.hide()

    def on_checkbox_toggled(self, state):
        if self.auto_lin_range.isChecked():
            self.time_min.hide()
            self.time_max.hide()
            self.time_step.hide()
            self.time_min_label.hide()
            self.time_max_label.hide()
            self.time_step_label.hide()
        else:
            self.time_min.show()
            self.time_max.show()
            self.time_step.show()
            self.time_min_label.show()
            self.time_max_label.show()
            self.time_step_label.show()
    
    def closeEvent(self, event):
        work_dir = os.environ.get('WORKING_DIR')
        temp_files = ["substrate_data.txt",
                       "time_data.txt",
                       "mutant.txt",
                       "name_data.txt",
                       "data_type.txt",
                       "blank_rep_data.txt",
                       "path_rep_data.txt",
                       "path_data.txt"
                       "sites.txt"]
        for file in temp_files:
            if os.path.exists(os.path.join(work_dir, file)):
                os.remove(os.path.join(work_dir, file))
            else:
                print("File already deleted")
        event.accept()

if __name__ == "__main__":
    try:
        app = QApplication(sys.argv)
        window = KineticAnalysisTool()
        window.show()
        sys.exit(app.exec_())
    except Exception as e:
        print("An error occurred:", e)

