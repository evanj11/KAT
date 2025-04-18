#!/usr/bin/env python
# coding: utf-8

# In[3]:


import customtkinter as ctk
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
import subprocess
from PIL import Image,ImageTk
import os


# In[4]:


def graph_hill_averages_winset():
    try:
        subprocess.run(['python', 'hill-averages_winset.py'], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running script: {e}")
    except FileNotFoundError:
        print("Script not found, Please check the path.")

def graph_hill_averages_noinset():
    try:
        subprocess.run(['python', 'hill-averages_noinset.py'], check=True)
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
        subprocess.run(['python', 'mm-averages.py'], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running script: {e}")
    except FileNotFoundError:
        print("Script not found, Please check the path.")

def graph_competitive():
    def comp_aver():
        try:
            subprocess.run(['python', 'comp-averages.py'], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running script: {e}")
        except FileNotFoundError:
            print("Script not found, Please check the path.")
    def comp_bf():
        try:
            subprocess.run(['python', 'comp-best_fit.py'], check=True)
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
        subprocess.run(['python', 'hill-averages_winset.py'], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running script: {e}")
    except FileNotFoundError:
        print("Script not found, Please check the path.")

def graph_noncompetitive():
    try:
        subprocess.run(['python', 'hill-averages_winset.py'], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running script: {e}")
    except FileNotFoundError:
        print("Script not found, Please check the path.")

def disp_name_entry():
    root = ctk.CTk()
    root.title("Name?")
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
    photo = ImageTk.PhotoImage(image)
    label = ctk.CTkImage(im, image=photo)
    label.image = photo 
    label.pack()

def disp_substrate_entry():
    root = ctk.CTk()
    root.title("Substrate?")
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
    root = ctk.CTk()
    root.title("Substrate?")
    time_min = ctk.CTkEntry(root)
    time_min_label = ctk.CTkLabel(root, text="Minimum Time to Read Data (min)")
    time_min.grid(row=1, column=1, pady=10)
    time_min_label.grid(row=2, column=1, pady=10)
    time_max = ctk.CTkEntry(root)
    time_max_label = ctk.CTkLabel(root, text="Maximum Time to Read Data (min)")
    time_max.grid(row=1, column=2, pady=10)
    time_max_label.grid(row=2, column=2, pady=10)
    step = ctk.CTkEntry(root)
    step_label = ctk.CTkLabel(root, text="Time-Step to Compute Slopes")
    step.grid(row=1, column=3, pady=10)
    step_label.grid(row=2, column=3, pady=10)
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
    submit_btn.grid(row=3, column=2, pady=10)


# In[41]:


class ColumnSelectorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Kinetic Analysis Tool (KAT)")
        self.columns = []
        for i in range(1,13):
            val = f"A{i}"
            self.columns.append(val)
        for i in range(1,13):
            val = f"B{i}"
            self.columns.append(val)
        self.selected_columns = []

        self.tree = ttk.Treeview(root, columns=self.columns, show="headings", selectmode="none")
        self.tree.grid(row=3, column=2, columnspan=10, pady=10)

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
        self.label.grid(row=4, column=5, pady=10)

    def toggle_column(self, index):
        if index in self.selected_columns:
            self.selected_columns.remove(index)
        else:
            self.selected_columns.append(index)
        self.label.config(text=f"Selected columns: {self.selected_columns}")
        with open("column_data.txt", "w") as file:
            file.writelines(f"{item}\n" for item in self.selected_columns)


# In[42]:


ctk.set_appearance_mode("System")
ctk.set_default_color_theme("blue")

window = ctk.CTk()
window.geometry("1900x700+100+200")
window.title("Kinetic Analysis Tool (KAT)")
my_font = ctk.CTkFont(family="Times New Roman", size=11)

#notebook = ttk.Notebook(window)

#tab1 = ttk.Frame(notebook)
#notebook.add(tab1, text="Hill Kinetics")

#tab2 = ttk.Frame(notebook)
#notebook.add(tab2, text="Michelis-Menten Kinetics")

#tab3 = ttk.Frame(notebook)
#notebook.add(tab3, text="Inhibition Kinetics")

tabview = ctk.CTkTabview(window)
tabview.grid(column=13, row=1, rowspan=3, columnspan=3, padx=20, pady=10)

tabview.add("Hill Kinetics")
tabview.add("Michelis-Menten Kinetics")
tabview.add("Inhibition Kinetics")

#kat = Image.open("/home/evanj/Downloads/ChatGPT Image Apr 16, 2025, 12_38_37 PM.png")
kat = Image.open("/home/evanj/Downloads/ChatGPT Image Apr 16, 2025, 04_48_48 PM.png")
resized_image = kat.resize((200,250))
photo_cat = ImageTk.PhotoImage(resized_image)
cat = ctk.CTkLabel(window, image=photo_cat)
cat.grid(row=6, column=14, rowspan=4, pady=10)

run_button1 = ctk.CTkButton(tabview.tab("Hill Kinetics"), text="Graph Averages Hill", command=plot_with_inset, corner_radius=12, font=my_font)
run_button1.grid(row=1, column=1, pady=10, padx=15)

run_button1a = ctk.CTkButton(tabview.tab("Hill Kinetics"), text="Graph Best Fit Hill", command=graph_hill_best_fit, corner_radius=12, font=my_font)
run_button1a.grid(row=1, column=2, pady=10, padx=15)

run_button2 = ctk.CTkButton(tabview.tab("Michelis-Menten Kinetics"), text="Graph Averages M-M", command=graph_mm_averages, corner_radius=12, font=my_font)
run_button2.grid(row=1, column=1, pady=10, padx=15)

run_button2a = ctk.CTkButton(tabview.tab("Michelis-Menten Kinetics"), text="Graph Best Fit M-M", command=graph_mm_best_fit, corner_radius=12, font=my_font)
run_button2a.grid(row=1, column=2, pady=10, padx=15)

run_button3 = ctk.CTkButton(tabview.tab("Inhibition Kinetics"), text="Competitive?", command=graph_competitive, corner_radius=12, font=my_font)
run_button3.grid(row=1, column=1, pady=10, padx=3)

run_button3a = ctk.CTkButton(tabview.tab("Inhibition Kinetics"), text="Uncompetitve?", command=graph_uncompetitive)
run_button3a.grid(row=1, column=2, pady=10, padx=3)

run_button3b = ctk.CTkButton(tabview.tab("Inhibition Kinetics"), text="Non-Competitive?", command=graph_noncompetitive)
run_button3b.grid(row=1, column=3, pady=10, padx=3)

graph_btn = ctk.CTkButton(window, text="Display Graph", command=disp_graph)
graph_btn.grid(row=5, column=13, pady=10)

name_btn = ctk.CTkButton(window, text="Name", command=disp_name_entry)
name_btn.grid(row=1, column=1, pady=10)
name = ctk.StringVar()
name_label = ctk.CTkLabel(window, textvariable=name)
name_label.grid(row=1, column=3, pady=10)

path_btn = ctk.CTkButton(window, text="File Path", command=open_file_explorer)
path_btn.grid(row=2, column=1, pady=10)

selected_path = ctk.StringVar()
path_label = ctk.CTkLabel(window, textvariable=selected_path)
path_label.grid(row=2, column=3, pady=10)

app = ColumnSelectorApp(window)

sub_btn = ctk.CTkButton(window, text="Substrate Information", command=disp_substrate_entry)
sub_btn.grid(row=6, column=1, pady=10)
peak = ctk.StringVar()
num = ctk.StringVar()
fac = ctk.StringVar()
sub_label1 = ctk.CTkLabel(window, textvariable=num)
sub_label1.grid(row=6, column=2, pady=10)
sub_label2 = ctk.CTkLabel(window, textvariable=fac)
sub_label2.grid(row=6, column=3, pady=10)
sub_label3 = ctk.CTkLabel(window, textvariable=peak)
sub_label3.grid(row=6, column=4, pady=10)
sub_labela = ctk.CTkLabel(window, text="Number of Substrates")
sub_labela.grid(row=5, column=2, pady=10)
sub_labelb = ctk.CTkLabel(window, text="Dilution Factor")
sub_labelb.grid(row=5, column=3, pady=10)
sub_labelc = ctk.CTkLabel(window, text="Maximum Substrate Conc.")
sub_labelc.grid(row=5, column=4, pady=10)

time_btn = ctk.CTkButton(window, text="Time Information", command=disp_time_entry)
time_btn.grid(row=8, column=1, pady=10)
min1 = ctk.StringVar()
max1 = ctk.StringVar()
step1 = ctk.StringVar()
sub_label1 = ctk.CTkLabel(window, textvariable=min1)
sub_label1.grid(row=8, column=2, pady=10)
sub_label2 = ctk.CTkLabel(window, textvariable=max1)
sub_label2.grid(row=8, column=3, pady=10)
sub_label3 = ctk.CTkLabel(window, textvariable=step1)
sub_label3.grid(row=8, column=4, pady=10)
sub_labela = ctk.CTkLabel(window, text="Minimum Time Value")
sub_labela.grid(row=7, column=2, pady=10)
sub_labelb = ctk.CTkLabel(window, text="Maximum Time Value")
sub_labelb.grid(row=7, column=3, pady=10)
sub_labelc = ctk.CTkLabel(window, text="Time-Step")
sub_labelc.grid(row=7, column=4, pady=10)



#notebook.grid(column=13, row=0, rowspan=3, columnspan=5, padx=20, pady=10)


# In[43]:


window.mainloop()


# In[35]:


os.remove("entry_data.txt")
os.remove("path_data.txt")
os.remove("column_data.txt")
os.remove("substrate_data.txt")


# In[ ]:




