from tkinter import *
import tkinter as tk
from tkinter import filedialog as fd
from treetracer import *

root = tk.Tk()
# canvas = tk.Canvas(root, width=600, height=600)
# canvas.pack()
root.geometry("900x850")

my_label = Label(root, text='Input FASTA and Newick file then choose analysis', font=("Helvetica", 19))
my_label.pack(pady = 40)

text_file_fasta = None
text_file_newick = None


def open_txt():
    global text_file_fasta
    text_file_fasta = fd.askopenfilename(initialdir="/Documents/Morton-Research/2021-research/", title="Open Text File", filetypes=(("Text Files", "*txt"),))
    text_file = open(text_file_fasta, 'r')
    stuff = text_file.read()
    my_text.insert(END, stuff)
    text_file.close()


def open_txt2():
    global text_file_newick
    text_file_newick = fd.askopenfilename(initialdir="/Documents/Morton-Research/2021-research/", title="Open Text File", filetypes=(("Text Files", "*txt"),))
    text_file = open(text_file_newick, 'r')
    stuff = text_file.read()
    my_text2.insert(END, stuff)
    text_file.close()


my_text = Text(root, width=40, height=5, font=("Helvetica", 13))
my_text.pack()

open_button = Button(root, text="Open FASTA File", command=open_txt)
open_button.pack()

my_text2 = Text(root, width=40, height=5, font=("Helvetica", 13))
my_text2.pack()

open_button = Button(root, text="Open Newick File", command=open_txt2)
open_button.pack()

def submit():
    # grab what is in my_box and pass into function
    my_text3.delete('1.0', END)
    if text_file_newick and text_file_fasta:
        tree_obj = TreeTracer(text_file_newick, text_file_fasta)
    else:
        tree_obj = TreeTracer('../iqtree_newick.txt', '../grass_rbcl_nodes_seq_fasta.txt')
    program = my_box.get()
    if program == 'n0':
        output = tree_obj.trace_tree_function(n0_context, branch_length=False)
    if program == 'n1':
        output = tree_obj.trace_tree_function(n1_context, branch_length=False)
    if program == 'n2':
        output = tree_obj.trace_tree_function(n2_context, branch_length=False)
    my_text3.insert(END, output)

my_box = Entry(root, fg="black", bg="lightblue", width=20)
my_box.pack(pady=30)

my_button = Button(root, text = 'RUN: Analysis Choice --> n0, n1, n2', command=submit)
my_button.pack()

my_text3 = Text(root, width=80, height=20, font=("Helvetica", 13))
my_text3.pack()

root.mainloop()



# label = tk.Label(text="Name")
# entry = tk.Entry()
# label.pack()
# entry.pack()
#
# name = entry.get()
# #
# # text_box = tk.Text()
# # text_box.pack()
# # words = text_box.get("1.0", tk.END)
#
# frame1 = tk.Frame(master=root, width=100, height=100, bg="red")
# frame1.pack(fill = tk.Y)
#