from tkinter import *
import tkinter as tk
from tkinter import filedialog as fd
from treetracer import *

root = tk.Tk()
root.geometry("900x750")
# bg=tk.PhotoImage(file='../background3.png')
# background_label = tk.Label( image=bg)
# background_label.place(x=0, y=0, relwidth=1, relheight=1)


my_label = Label(root, text='Input FASTA and Newick file then choose analysis', font=("Helvetica", 19),
                 background='#BABAD3')
my_label.pack(pady=10)

text_file_fasta = None
text_file_newick = None


def open_txt():
    global text_file_fasta
    text_file_fasta = fd.askopenfilename(initialdir="/Documents/Morton-Research/2021-research/", title="Open Text File",
                                         filetypes=(("Text Files", "*txt"),))
    text_file = open(text_file_fasta, 'r')
    stuff = text_file.read()
    my_text.insert(END, stuff)
    text_file.close()


def open_txt2():
    global text_file_newick
    text_file_newick = fd.askopenfilename(initialdir="/Documents/Morton-Research/2021-research/",
                                          title="Open Text File", filetypes=(("Text Files", "*txt"),))
    text_file = open(text_file_newick, 'r')
    stuff = text_file.read()
    my_text2.insert(END, stuff)
    text_file.close()


my_text = Text(root, width=60, height=10, font=("Helvetica", 15), bd=5, relief=GROOVE)
my_text.pack()

open_button = Button(root, text="Open FASTA File", command=open_txt, activeforeground='blue')
open_button.pack()

my_text2 = Text(root, width=60, height=5, font=("Helvetica", 15), bd=5, relief=GROOVE)
my_text2.pack()

open_button = Button(root, text="Open Newick File", command=open_txt2, activeforeground='blue')
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
        o2 = tree_obj.print_cumulative_matrices()
    if program == 'n1':
        output = tree_obj.trace_tree_function(n1_context, branch_length=False)
        o2 = tree_obj.print_cumulative_matrices()
    if program == 'n2':
        output = tree_obj.trace_tree_function(n2_context, branch_length=False)
        o2 = tree_obj.print_cumulative_matrices()
    if program == 'n0 fourfold':
        output = tree_obj.trace_tree_function(fourfold_n0_context, branch_length=False)
        o2 = tree_obj.print_cumulative_matrices()
    if program == 'n1 fourfold':
        output = tree_obj.trace_tree_function(fourfold_n1_context, branch_length=False)
        o2 = tree_obj.print_cumulative_matrices()
    if program == 'n2 fourfold':
        output = tree_obj.trace_tree_function(fourfold_n2_context, branch_length=False)
        o2 = tree_obj.print_cumulative_matrices()
    my_text3.insert(END, output)
    my_text3.insert(END, o2)


label = Label(root, text='Enter analysis option below:\n Analysis Choice --> n0, n1, n2, n0 fourfold, '
                         'n1 fourfold, n2 fourfold', relief=RAISED)
label.pack(pady=30)

my_box = Entry(root, fg="black", bg="lightblue", width=20, relief=GROOVE)
my_box.pack()

my_button = Button(root, text='RUN', command=submit, activeforeground='blue')
my_button.pack(pady=10)

my_text3 = Text(root, width=50, height=10, font=("Helvetica", 15), bd=2, relief=GROOVE)
my_text3.pack(pady=10)

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
