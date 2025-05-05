#!/usr/bin/env python
# coding: utf-8

# In[1]:


import customtkinter as ctk
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
import subprocess
from PIL import Image,ImageTk
import os


# In[2]:


def graph_hill_averages_winset():
    try:
        subprocess.run(['python', '../scripts/hill-averages_winset.py'], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running script: {e}")
    except FileNotFoundError:
        print("Script not found, Please check the path.")

def graph_hill_averages_noinset():
    try:
        subprocess.run(['python', '../scripts/hill-averages_noinset.py'], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running script: {e}")
    except FileNotFoundError:
        print("Script not found, Please check the path.")

def plot_with_inset():
    root = ctk.CTk()
    root.title("Inset Plot?")
    yes = ctk.CTkButton(root, text='Yes', command=graph_hill_averages_winset)
    yes.grid(row=2, column=1, pady=10, padx=20)
    no = ctk.CTkButton(root, text='No', command=graph_hill_averages_noinset)
    no.grid(row=2, column=2, pady=10, padx=20)
    n_label = ctk.CTkLabel(root, text="Do you want to include an inset plot?")
    n_label.grid(row=1, column=1, columnspan=2, pady=10)
    def close():
        root.destroy()
    close_button = ctk.CTkButton(root, text='Close', command=close)
    close_button.grid(row=3, column=1, columnspan=2, pady=10)

def graph_mm_averages():
    try:
        subprocess.run(['python', '../scripts/mm-averages.py'], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running script: {e}")
    except FileNotFoundError:
        print("Script not found, Please check the path.")

def graph_competitive():
    def comp_aver():
        try:
            subprocess.run(['python', '../scripts/comp-averages.py'], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running script: {e}")
        except FileNotFoundError:
            print("Script not found, Please check the path.")
    def comp_bf():
        try:
            subprocess.run(['python', '../scripts/comp-best_fit.py'], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running script: {e}")
        except FileNotFoundError:
            print("Script not found, Please check the path.")
    root = ctk.CTk()
    root.title("Plot Competitive Inhibition")
    average = ctk.CTkButton(root, text='Graph Averages?', command=comp_aver)
    average.grid(row=2, column=1, pady=10, padx=20)
    bf = ctk.CTkButton(root, text='Graph Best Fit?', command=comp_bf)
    bf.grid(row=2, column=2, pady=10, padx=20)
    n_label = ctk.CTkLabel(root, text="How do you want to compute velocities?")
    n_label.grid(row=1, column=1, columnspan=2, pady=10)
    def close():
        root.destroy()
    close_button = ctk.CTkButton(root, text='Close', command=close)
    close_button.grid(row=3, column=1, columnspan=2, pady=10)    

def graph_uncompetitive():
    try:
        subprocess.run(['python', '../scritps/hill-averages_winset.py'], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running script: {e}")
    except FileNotFoundError:
        print("Script not found, Please check the path.")

def graph_noncompetitive():
    try:
        subprocess.run(['python', '../scripts/hill-averages_winset.py'], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running script: {e}")
    except FileNotFoundError:
        print("Script not found, Please check the path.")

def disp_name_entry():
    root = ctk.CTkToplevel()
    root.title("Name?")
    root.grab_set() 
    entry = ctk.CTkEntry(root)
    entry.grid(row=2, column=1, pady=10)
    def submit():
        data = entry.get()
        name.set(data)
        with open("entry_data.txt", "w") as f:
            f.write(data)
        print("Entry saved.")
        root.destroy()
    submit_btn = ctk.CTkButton(root, text="Submit", command=submit)
    submit_btn.grid(row=3, column=1, pady=10)
    n_label = ctk.CTkLabel(root, text="Name of Output File")
    n_label.grid(row=1, column=1, pady=10)

def open_file_explorer():
    file_path = filedialog.askopenfilename()
    if file_path:
        selected_path.set(file_path)
        with open("path_data.txt", "w") as f:
            f.write(file_path)

def disp_graph():
    with open("entry_data.txt", 'r') as file:
        file = file.readlines()
    filename = file[0]
    im = tk.Toplevel()
    image = Image.open(f"{filename}.png")
    photo = ctk.CTkImage(light_image=image, size=(1000,800))
    label = ctk.CTkLabel(im, image=photo, text="")
    label.image = photo 
    label.pack()

def disp_substrate_entry():
    root = ctk.CTkToplevel()
    root.title("Substrate?")
    root.grab_set()
    dil_num = ctk.CTkEntry(root)
    dil_num_label = ctk.CTkLabel(root, text="Number of Substrate Conc.")
    dil_num.grid(row=1, column=1, pady=10)
    dil_num_label.grid(row=2, column=1, pady=10)
    dil = ctk.CTkEntry(root)
    dil_label = ctk.CTkLabel(root, text="Dilution Factor")
    dil.grid(row=1, column=2, pady=10)
    dil_label.grid(row=2, column=2, pady=10)
    dil_peak = ctk.CTkEntry(root)
    dil_peak_label = ctk.CTkLabel(root, text="Highest Substrate Conc.")
    dil_peak.grid(row=1, column=3, pady=10)
    dil_peak_label.grid(row=2, column=3, pady=10)
    def submit():
        data1 = dil_num.get()
        num.set(data1)
        data2 = dil.get()
        fac.set(data2)
        data3 = dil_peak.get()
        peak.set(data3)
        with open("substrate_data.txt", "w") as f:
            f.write(f"{data1}\n")
            f.write(f"{data2}\n")
            f.write(f"{data3}\n")
        print("Entry saved.")
        root.destroy()
    submit_btn = ctk.CTkButton(root, text="Submit", command=submit)
    submit_btn.grid(row=3, column=2, pady=10)

def graph_hill_best_fit():
    try:
        subprocess.run(['python', 'hill-best_fit.py'], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running script: {e}")
    except FileNotFoundError:
        print("Script not found, Please check the path.")

def graph_mm_best_fit():
    try:
        subprocess.run(['python', 'mm-best_fit.py'], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running script: {e}")
    except FileNotFoundError:
        print("Script not found, Please check the path.")


def disp_time_entry():
    root = ctk.CTkToplevel()
    root.title("Time-Course?")
    root.grab_set()
    time_min = ctk.CTkEntry(root)
    time_min_label = ctk.CTkLabel(root, text="Minimum Time to Read Data (min)")
    time_min.grid(row=1, column=1, pady=10)
    time_min_label.grid(row=2, column=1, pady=10, padx=10)
    time_max = ctk.CTkEntry(root)
    time_max_label = ctk.CTkLabel(root, text="Maximum Time to Read Data (min)")
    time_max.grid(row=1, column=2, pady=10)
    time_max_label.grid(row=2, column=2, pady=10, padx=10)
    step = ctk.CTkEntry(root)
    step_label = ctk.CTkLabel(root, text="Time-Step to Compute Slopes")
    step.grid(row=1, column=3, pady=10)
    step_label.grid(row=2, column=3, pady=10, padx=10)
    label = ctk.CTkLabel(root, text="*Difference in Max and Min time must be divisible by Time-Step*")
    label.grid(row=3, column=1, columnspan=3, pady=10)
    label.configure(text_color="red")
    def submit():
        data1 = time_min.get()
        min1.set(data1)
        data2 = time_max.get()
        max1.set(data2)
        data3 = step.get()
        step1.set(data3)
        with open("time_data.txt", "w") as f:
            f.write(f"{data1}\n")
            f.write(f"{data2}\n")
            f.write(f"{data3}\n")
        print("Entry saved.")
        root.destroy()
    submit_btn = ctk.CTkButton(root, text="Submit", command=submit)
    submit_btn.grid(row=4, column=2, pady=10)

def abs_flo():
    root = ctk.CTkToplevel()
    root.title("Absorbance or Fluorescence?")
    root.grab_set()
    def absor():
        dat_type.set("Absorbance \n")
        root = ctk.CTkToplevel()
        root.title("Name?")
        root.grab_set()
        entry = ctk.CTkEntry(root)
        entry.grid(row=2, column=1, pady=10)
        def submit():
            data = entry.get()
            molabs.set(data)
            with open("data_type.txt", "w") as f:
                f.write("absorbance\n")
                f.write(f"{data}\n")
            root.destroy()
        submit_btn = ctk.CTkButton(root, text="Submit", command=submit)
        submit_btn.grid(row=3, column=1, pady=10)
        n_label = ctk.CTkLabel(root, text="Molar Absorptivity")
        n_label.grid(row=1, column=1, pady=10)
    def fluo():
        dat_type.set("Fluorescence")
        with open("data_type.txt", "w") as f:
            f.write("fluorescence")
        root.destroy()
    def ok():
            root.destroy()
    absorb = ctk.CTkButton(root, text='Absorbance', command=absor)
    absorb.grid(row=2, column=1, pady=10, padx=20)
    fluor = ctk.CTkButton(root, text='Fluorescence', command=fluo)
    fluor.grid(row=2, column=2, pady=10, padx=20)
    ok_btn = ctk.CTkButton(root, text="Ok", command=ok)
    ok_btn.grid(row=3, column=1, pady=10, columnspan=2)
    n_label = ctk.CTkLabel(root, text="Are you using absorbance or fluorescence data?")
    n_label.grid(row=1, column=1, columnspan=2, pady=10)


# In[3]:


class ColumnSelectorApp:
    def __init__(self, root):
        self.root = root
        #self.root.title("Kinetic Analysis Tool (KAT)")
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)
        self.columns = []
        for i in range(1,13):
            val = f"A{i}"
            self.columns.append(val)
        for i in range(1,13):
            val = f"B{i}"
            self.columns.append(val)
        self.selected_columns = []

        self.tree = ttk.Treeview(root, height=6, columns=self.columns, show="headings", selectmode="none")
        self.tree.grid(row=3, column=1, columnspan=10, rowspan=2, pady=10)

        # Create styles for selected and normal column headers
        self.style = ttk.Style()
        self.style.configure("Treeview.Heading", font=('Segoe UI', 10))
        self.style.configure("Selected.Treeview.Heading", background="#d0e0ff")

        # Setup headers with command for click behavior
        for idx, col in enumerate(self.columns):
            self.tree.heading(col, text=col, command=lambda i=idx: self.toggle_column(i))
            self.tree.column(col, width=50)

        # Label to show selection
        self.label = ctk.CTkLabel(root, text="Selected columns: []")
        self.label.grid(row=5, column=5, pady=10)

    def toggle_column(self, index):
        if index in self.selected_columns:
            self.selected_columns.remove(index)
        else:
            self.selected_columns.append(index)
        self.label.configure(text=f"Selected columns: {self.selected_columns}", text_color="green")
        with open("column_data.txt", "w") as file:
            file.writelines(f"{item}\n" for item in self.selected_columns)


# In[4]:


ctk.set_appearance_mode("System")
ctk.set_default_color_theme("dark-blue")

#screen = ctk.CTk()
#window = ctk.CTkFrame(master=screen)
#window.grid(row=0, column=0)

window = ctk.CTk()
window.state("normal")
window.title("Kinetic Analysis Tool (KAT)")
window.grid_rowconfigure(0, weight=1)
window.grid_columnconfigure(0, weight=1)
my_font = ctk.CTkFont(family="Times New Roman", size=18)

tabview = ctk.CTkTabview(window)
tabview.add("Hill Kinetics")
tabview.add("Michaelis-Menten Kinetics")
tabview.add("Inhibition Kinetics")
tabview.grid(column=13, row=1, rowspan=3, padx=20, pady=10)

kat = Image.open("../.imgs/ChatGPT Image Apr 18, 2025, 12_39_18 PM.png")
photo_cat = ctk.CTkImage(light_image=kat, size=(325,425))
cat = ctk.CTkLabel(window, image=photo_cat, text="")
cat.grid(row=11, column=13, rowspan=4, columnspan=5, pady=10)
kat1 = Image.open("../.imgs/ChatGPT Image Apr 18, 2025, 12_39_18 PM.png")
photo_cat1 = ctk.CTkImage(light_image=kat, size=(325,425))
cat1 = ctk.CTkLabel(window, image=photo_cat, text="")
cat1.grid(row=11, column=0, rowspan=1, columnspan=2, pady=10)

gui_label = ctk.CTkLabel(window, anchor=ctk.W, font=my_font, text="Welcome to KAT, the Kinetic Analysis Toolkit! \n To generate a kinetic curve following either Hill or Michelis-Menten kinetics, follow these steps: \n 1. Create a .csv file with values for fluorescence or absorbance that follows [Time] [Temp] [Data]. \n 2. Give the output graph file a name. \n 3. Specify a path to the data in a .csv format. \n 4. Specify the kind of data you have collected (absorbance vs. fluorescence). \n 5. Input the information about the substrate concentrations used (must be serially diluted). \n 6. Input the information about the time-course to sample. \n 7. Select the columns to read data from. \n 8. Click either *Graph Averages* or *Graph Best-Fit* for Hill or Michelis-Menten Kinetics. \n 9. Click *Display Graph* to open the graph output.")
gui_label.grid(row=11, column=2, rowspan=1, columnspan=10, pady=10)

run_button1 = ctk.CTkButton(tabview.tab("Hill Kinetics"), text="Graph Averages Hill", command=plot_with_inset)
run_button1.grid(row=1, column=1, columnspan=2, pady=15, padx=15)

run_button1a = ctk.CTkButton(tabview.tab("Hill Kinetics"), text="Graph Best Fit Hill", command=graph_hill_best_fit)
run_button1a.grid(row=2, columnspan=2, column=1, pady=15, padx=15)

run_button2 = ctk.CTkButton(tabview.tab("Michaelis-Menten Kinetics"), text="Graph Averages M-M", command=graph_mm_averages)
run_button2.grid(row=1, column=1, columnspan=1, pady=15, padx=15)

run_button2a = ctk.CTkButton(tabview.tab("Michaelis-Menten Kinetics"), text="Graph Best Fit M-M", command=graph_mm_best_fit)
run_button2a.grid(row=2, column=1, columnspan=1, pady=15, padx=15)

mm = Image.open("../.imgs/mm.png")
photo_mm = ctk.CTkImage(light_image=mm, size=(155,55))
mmeq = ctk.CTkLabel(tabview.tab("Michaelis-Menten Kinetics"), image=photo_mm, text="")
mmeq.grid(row=1, column=2, rowspan=4, columnspan=2, pady=15, padx=15)

run_button3 = ctk.CTkButton(tabview.tab("Inhibition Kinetics"), text="Competitive?", command=graph_competitive)
run_button3.grid(row=1, column=1, columnspan=2, pady=15, padx=15)

run_button3a = ctk.CTkButton(tabview.tab("Inhibition Kinetics"), text="Uncompetitve?", command=graph_uncompetitive)
run_button3a.grid(row=2, column=1, columnspan=2, pady=15, padx=15)

run_button3b = ctk.CTkButton(tabview.tab("Inhibition Kinetics"), text="Non-Competitive?", command=graph_noncompetitive)
run_button3b.grid(row=3, column=1, columnspan=2, pady=15, padx=15)

graph_btn = ctk.CTkButton(window, text="Display Graph", command=disp_graph)
graph_btn.grid(row=4, column=13, pady=10)

name_btn = ctk.CTkButton(window, text="Name", command=disp_name_entry)
name_btn.grid(row=1, column=0, pady=10)
name = ctk.StringVar()
name_label = ctk.CTkLabel(window, textvariable=name)
name_label.grid(row=1, column=2, pady=10)
name_label.configure(text_color="green")

path_btn = ctk.CTkButton(window, text="File Path", command=open_file_explorer)
path_btn.grid(row=2, column=0, pady=10)

selected_path = ctk.StringVar()
path_label = ctk.CTkLabel(window, textvariable=selected_path)
path_label.grid(row=2, column=1, columnspan=3, pady=10)
path_label.configure(text_color="green")

app = ColumnSelectorApp(window)

dat_btn = ctk.CTkButton(window, text="Select Data Type", command=abs_flo)
dat_btn.grid(row=6, column=0, pady=10)
dat_type = ctk.StringVar()
dat_label = ctk.CTkLabel(window, textvariable=dat_type)
dat_label.grid(row=6, column=1, pady=10)
dat_label.configure(text_color="green")
dat_labela = ctk.CTkLabel(window, text="Data Type")
dat_labela.grid(row=5, column=1, pady=10)
molabs = ctk.StringVar()
molabs_label = ctk.CTkLabel(window, textvariable=molabs)
molabs_label.grid(row=6, column=2, pady=10)
molabs_label.configure(text_color="green")
molabs_labela = ctk.CTkLabel(window, text="Molar Absorptivity")
molabs_labela.grid(row=5, column=2, pady=10)

sub_btn = ctk.CTkButton(window, text="Substrate Information", command=disp_substrate_entry)
sub_btn.grid(row=8, column=0, pady=10, padx=10)
peak = ctk.StringVar()
num = ctk.StringVar()
fac = ctk.StringVar()
sub_label1 = ctk.CTkLabel(window, textvariable=num)
sub_label1.grid(row=8, column=1, pady=10)
sub_label1.configure(text_color="green")
sub_label2 = ctk.CTkLabel(window, textvariable=fac)
sub_label2.grid(row=8, column=2, pady=10)
sub_label2.configure(text_color="green")
sub_label3 = ctk.CTkLabel(window, textvariable=peak)
sub_label3.grid(row=8, column=3, pady=10)
sub_label3.configure(text_color="green")
sub_labela = ctk.CTkLabel(window, text="Number of Substrates")
sub_labela.grid(row=7, column=1, pady=10)
sub_labelb = ctk.CTkLabel(window, text="Dilution Factor")
sub_labelb.grid(row=7, column=2, pady=10)
sub_labelc = ctk.CTkLabel(window, text="Maximum Substrate Conc.")
sub_labelc.grid(row=7, column=3, pady=10)

time_btn = ctk.CTkButton(window, text="Time Information", command=disp_time_entry)
time_btn.grid(row=10, column=0, pady=10, padx=5)
min1 = ctk.StringVar()
max1 = ctk.StringVar()
step1 = ctk.StringVar()
time_label1 = ctk.CTkLabel(window, textvariable=min1)
time_label1.grid(row=10, column=1, pady=10)
time_label1.configure(text_color="green")
time_label2 = ctk.CTkLabel(window, textvariable=max1)
time_label2.grid(row=10, column=2, pady=10)
time_label2.configure(text_color="green")
time_label3 = ctk.CTkLabel(window, textvariable=step1)
time_label3.grid(row=10, column=3, pady=10)
time_label3.configure(text_color="green")
time_labela = ctk.CTkLabel(window, text="Minimum Time Value")
time_labela.grid(row=9, column=1, pady=10)
time_labelb = ctk.CTkLabel(window, text="Maximum Time Value")
time_labelb.grid(row=9, column=2, pady=10)
time_labelc = ctk.CTkLabel(window, text="Time-Step")
time_labelc.grid(row=9, column=3, pady=10)

window.mainloop()


# In[72]:


os.remove("entry_data.txt")
os.remove("path_data.txt")
os.remove("column_data.txt")
os.remove("substrate_data.txt")
os.remove("data_type.txt")
os.remove("time_data.txt")


# In[ ]:




