from PyQt5.QtWidgets import QMainWindow, QApplication, QLabel, QMessageBox, QCheckBox, QPushButton, QFileDialog, QPlainTextEdit
from PyQt5 import uic
import sys
from treetracer import TreeTracer
from treetracer import *


class TreeTracerUI(QMainWindow):
    def __init__(self):
        super(TreeTracerUI, self).__init__()

        # load the .ui file
        uic.loadUi("TreeTracerUI_test.ui", self)

        # define widgets
        self.newickfile_label = self.findChild(QLabel, "newickfile_label")
        self.newick_button = self.findChild(QPushButton, "newick_button")
        self.fastafile_label = self.findChild(QLabel, "fastafile_label")
        self.fasta_button = self.findChild(QPushButton, "fasta_button")
        self.newick_text = self.findChild(QPlainTextEdit, "newick_text")
        self.fasta_text = self.findChild(QPlainTextEdit, "fasta_text")
        self.n0 = self.findChild(QCheckBox, "n0")
        self.analysis_annotation = self.findChild(QLabel, "analysis_annotation")
        self.n_annotation = self.findChild(QLabel, "n_annotation")
        self.fourfold_anno = self.findChild(QLabel, "fourfold_anno")
        self.n1 = self.findChild(QCheckBox, "n1")
        self.n2 = self.findChild(QCheckBox, "n2")
        self.n0_fourfold = self.findChild(QCheckBox, "n0_fourfold")
        self.n1_fourfold = self.findChild(QCheckBox,"n1_fourfold")
        self.n2_fourfold = self.findChild(QCheckBox,"n2_fourfold")

        self.save_anno = self.findChild(QLabel, "save_anno")
        self.save_matrices = self.findChild(QCheckBox, "save_matrices")
        self.print_all_results = self.findChild(QCheckBox, "print_all_results")
        self.save_site_df = self.findChild(QCheckBox,"save_site_df")
        self.save_graphs = self.findChild(QCheckBox,"save_graphs")

        self.sitebysite_analysis = self.findChild(QCheckBox,"sitebysite_analysis")

        self.press_to_run_button = self.findChild(QPushButton, "press_to_run_button")
        self.output_label = self.findChild(QPlainTextEdit, "output_label")
        self.outgroup_entry = self.findChild(QPlainTextEdit, "output_label")
        self.push_matrix_info = self.findChild(QPushButton, "push_matrix_info")
        self.push_sitexsite_info = self.findChild(QPushButton, "push_sitexsite_info")
        self.push_saveresults_info = self.findChild(QPushButton, "push_saveresults_info")


        # variables for running TreeTracer
        # newick file you open
        self.newick_path = ""
        # fasta file you open
        self.fasta_path = ""
        # ADD IN: if these paths are still "" then take what is in the box,
        # if nothing in box then return error

        # open files for self.newick_button and self.fasta_button
        self.newick_button.clicked.connect(self.open_dialog_box_newick)
        self.fasta_button.clicked.connect(self.open_dialog_box_fasta)

        # click buttons to learn more about program to run
        self.push_matrix_info.clicked.connect(lambda: self.clicked_more_info("matrix"))
        self.push_sitexsite_info.clicked.connect(lambda: self.clicked_more_info("site"))
        self.push_saveresults_info.clicked.connect(lambda: self.clicked_more_info("saveresults"))


        # Execute when press_to_run_button is clicked
        self.press_to_run_button.clicked.connect(self.execute)

        # show the app
        self.show()

    def execute(self):
        matrix_output = []
        # check that there is a path to the sequence and newick files
        if len(self.newick_path) > 2 and len(self.fasta_path) > 2:
            if self.outgroup_entry:
                outgroups_input = self.outgroup_entry.toPlainText()
                print("OUTGROUPS:", outgroups_input)
            self.tree_obj = TreeTracer(self.newick_path, self.fasta_path, outgroups=outgroups_input)
            if self.n2_fourfold.isChecked():
                self.tree_obj.trace_tree_function(fourfold_n2_context, branch_length=False)
                self.tree_obj.print_cumulative_matrices()
            if self.n1_fourfold.isChecked():
                self.tree_obj.trace_tree_function(fourfold_n1_context, branch_length=False)
                self.tree_obj.print_cumulative_matrices()
            if self.n0_fourfold.isChecked():
                self.tree_obj.trace_tree_function(fourfold_n0_context, branch_length=False)
                self.tree_obj.print_cumulative_matrices()
            if self.n0.isChecked():
                self.tree_obj.trace_tree_function(n0_context, branch_length=False)
                self.print_matrix_results()
            if self.n1.isChecked():
                self.tree_obj.trace_tree_function(n1_context, branch_length=False)
                self.print_matrix_results()
            if self.n2.isChecked():
                self.tree_obj.trace_tree_function(n2_context, branch_length=False)
                self.tree_obj.print_cumulative_matrices()
            if self.sitebysite_analysis.isChecked():
                self.tree_obj.site_trace_tree_function()
                self.tree_obj.site_change_analysis(to_csv=False, show_graphs=False, save_graphs=False, run_stats=False)
        #self.output_label.setPlainText("done")
        #self.output_label.appendPlainText("added")
        #self.save_site_df.setChecked(False)

    def print_matrix_results(self):
        if self.print_all_results.isChecked():
            for m in self.tree_obj.print_cumulative_matrices():
                self.output_label.appendPlainText(str(m))
            self.output_label.appendPlainText("--------\n")
            #self.output_label.appendPlainText(self.tree_obj.print_cumulative_matrices())


    def open_dialog_box_newick(self):
        # open a dialog box
        filename = QFileDialog.getOpenFileName()
        print(filename[0])
        self.newick_path = filename[0]
        # open file
        if len(self.newick_path) < 1 or self.newick_path[-3:] != 'txt':
            self.newick_text.setPlainText('NEWICK FILE .TXT NOT SELECTED')
            return
        with open(self.newick_path, "r") as f:
            self.newick_text.setPlainText(f.read())

    def open_dialog_box_fasta(self):
        # open a dialog box
        filename = QFileDialog.getOpenFileName()
        print(filename[0])
        self.fasta_path = filename[0]
        if len(self.fasta_path) < 1 or self.fasta_path[-3:] != 'txt':
            self.fasta_text.setPlainText('FASTA SEQUENCES FILE .TXT NOT SELECTED')
            return
        # open file
        with open(self.fasta_path, "r") as f:
            self.fasta_text.setPlainText(f.read())

    def clicked_more_info(self, type_info):
        if type_info == "matrix":
            msg = QMessageBox()
            msg.setWindowTitle("Info regarding Matrix Analysis")
            msg.setText("n = number of surrounding bases for context dependent analysis\n"
                        "Fourfold = only counting mutation rate for 4-fold degenerate sites")
            x = msg.exec_()
        elif type_info == "site":
            msg = QMessageBox()
            msg.setWindowTitle("Info regarding Site by Site Analysis")
            msg.setText("--Site by Site Analysis--")
            x = msg.exec_()
        elif type_info == "saveresults":
            msg = QMessageBox()
            msg.setWindowTitle("Info regarding Saving and Printing Analysis Output")
            msg.setText("--saveresults--")
            x = msg.exec_()

# Initialize the app
app = QApplication(sys.argv)
UIWindow = TreeTracerUI()
app.exec_()
