# -*- coding: utf-8 -*-
# Clean UTF-8 safe GUI for InSilicoPCR

import tkinter as tk
from tkinter import filedialog, messagebox
import subprocess
import os

def browse_fasta_dir():
    path = filedialog.askdirectory()
    fasta_dir_var.set(path)

def browse_primers():
    path = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv")])
    primers_var.set(path)

def browse_output_dir():
    path = filedialog.askdirectory()
    output_dir_var.set(path)

def run_insilico_pcr():
    fasta_dir = fasta_dir_var.get()
    primers = primers_var.get()
    out_dir = output_dir_var.get()
    max_mm = max_mm_var.get()
    min_len = min_len_var.get()
    max_len = max_len_var.get()
    workers = workers_var.get()

    if not fasta_dir or not primers or not out_dir:
        messagebox.showerror("Error", "Please select all required paths.")
        return

    cmd = [
        "python",
        "insilico_pcr_parallel_multifasta.py",
        "--fasta_dir", fasta_dir,
        "--primers", primers,
        "--out_dir", out_dir,
        "--max_mismatch", str(max_mm),
        "--min_len", str(min_len),
        "--max_len", str(max_len),
        "--workers", str(workers)
    ]

    try:
        subprocess.run(cmd, check=True)
        messagebox.showinfo("Success", "InSilico PCR completed successfully.")
    except subprocess.CalledProcessError:
        messagebox.showerror("Error", "The tool encountered an error. Check parameters and try again.")

# GUI window setup
root = tk.Tk()
root.title("InSilicoPCR GUI - Dr. Vishal Mevada, NFSU")
root.geometry("600x450")

# Tkinter Variables
fasta_dir_var = tk.StringVar()
primers_var = tk.StringVar()
output_dir_var = tk.StringVar()
max_mm_var = tk.IntVar(value=2)
min_len_var = tk.IntVar(value=50)
max_len_var = tk.IntVar(value=2000)
workers_var = tk.IntVar(value=0)

# GUI Layout
tk.Label(root, text="Genome FASTA Directory:").pack()
tk.Entry(root, textvariable=fasta_dir_var, width=60).pack()
tk.Button(root, text="Browse", command=browse_fasta_dir).pack()

tk.Label(root, text="Primer CSV File:").pack()
tk.Entry(root, textvariable=primers_var, width=60).pack()
tk.Button(root, text="Browse", command=browse_primers).pack()

tk.Label(root, text="Output Directory:").pack()
tk.Entry(root, textvariable=output_dir_var, width=60).pack()
tk.Button(root, text="Browse", command=browse_output_dir).pack()

tk.Label(root, text="Max Mismatches:").pack()
tk.Entry(root, textvariable=max_mm_var).pack()

tk.Label(root, text="Min Amplicon Length:").pack()
tk.Entry(root, textvariable=min_len_var).pack()

tk.Label(root, text="Max Amplicon Length:").pack()
tk.Entry(root, textvariable=max_len_var).pack()

tk.Label(root, text="Workers (0 = Auto):").pack()
tk.Entry(root, textvariable=workers_var).pack()

tk.Button(root, text="Run InSilico PCR", command=run_insilico_pcr, bg="green", fg="white").pack(pady=10)

root.mainloop()

