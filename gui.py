from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QLineEdit, \
    QFileDialog, QProgressBar, QCheckBox, QTabWidget, QTextEdit
from processor import get_meta_data, get_exp_fps, get_highest_rep, get_exp_dupes, get_ratio, calc_chi
from main import bitreps_measure, dir_setup, RESULTS_DIR
from pathlib import Path
import sys
import os


def set_t2_stats_edit(exp="", obs="", chi="", efp="", ed="", oh="", ra="", mr=""):
    """
    Template string for automated statistical analysis display
    :param exp: Expected distribution
    :param obs: Observed distribution
    :param chi: Chi-square result
    :param efp: Expected number of false positives
    :param ed: Expected number of genuine duplicates
    :param oh: Observed number of repetitions
    :param ra: Num. expected hits / Num. observed hits
    :param mr: Maximally repeating block within RNG output
    :return: Formatted string of given parameters
    """
    return "Expected distribution: %s\n" \
           "Observed distribution: %s\n" \
           "Chi-square: %s\n\n" \
           "Expected false positives: %s\n" \
           "Expected duplicates: %s\n" \
           "Observed hits: %s\n" \
           "Ratio: %s\n\n" \
           "Maximum repetition: %s\n" % (exp, obs, chi, efp, ed, oh, ra, mr)


class BitReps(QWidget):
    def __init__(self):
        # Window setup
        super().__init__()
        self.width, self.height = 800, 200
        self.setMinimumSize(self.width, self.height)
        self.setWindowTitle("BitReps - University of Kent 2022 - Quantum Hub 2")

        # Main layout
        self.main_layout = QVBoxLayout()
        self.setLayout(self.main_layout)

        # Tab setup
        self.tabwidget = QTabWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()

        self.tab1_layout = QVBoxLayout()
        self.tab2_layout = QVBoxLayout()

        self.tab1.setLayout(self.tab1_layout)
        self.tab2.setLayout(self.tab2_layout)

        # ##### Tab 1 Contents ##### #
        # Input select button
        self.t1_input_btn = QPushButton("Select Input Data")
        self.t1_input_btn.clicked.connect(self.sel_t1_file)

        # Input select label and display
        self.t1_sub_1 = QHBoxLayout()
        self.t1_input_lab = QLabel("Chosen File:")
        self.t1_input_edit = QLineEdit()
        self.t1_input_edit.setPlaceholderText("Select File...")
        self.t1_input_edit.setReadOnly(True)
        self.t1_sub_1.addWidget(self.t1_input_lab)
        self.t1_sub_1.addWidget(self.t1_input_edit)

        # Blocksize label and input
        self.t1_sub_2 = QHBoxLayout()
        self.t1_size_lab = QLabel("Blocksize:")
        self.t1_size_edit = QLineEdit()
        self.t1_size_edit.setPlaceholderText("In bits... (16, 32, 64, 128, 256, 512)")
        self.t1_sub_2.addWidget(self.t1_size_lab)
        self.t1_sub_2.addWidget(self.t1_size_edit)

        # Sliding window label and checkbox
        self.t1_sub_3 = QHBoxLayout()
        self.t1_slide_lab = QLabel("Sliding window?:")
        self.t1_slide_chk = QCheckBox()
        self.t1_slide_chk.setEnabled(False)
        self.t1_sub_3.addWidget(self.t1_slide_lab)
        self.t1_sub_3.addWidget(self.t1_slide_chk)

        # Error rate label and input
        self.t1_sub_4 = QHBoxLayout()
        self.t1_err_lab = QLabel("Error rate:")
        self.t1_err_edit = QLineEdit()
        self.t1_err_edit.setPlaceholderText("Between 0 and 1...")
        self.t1_sub_4.addWidget(self.t1_err_lab)
        self.t1_sub_4.addWidget(self.t1_err_edit)

        # Reset and run buttons
        self.t1_sub_5 = QHBoxLayout()
        self.t1_rst_btn = QPushButton("Reset")
        self.t1_rst_btn.clicked.connect(self.t1_reset)
        self.t1_run_btn = QPushButton("Run")
        self.t1_run_btn.clicked.connect(self.measure)
        self.t1_sub_5.addWidget(self.t1_rst_btn)
        self.t1_sub_5.addWidget(self.t1_run_btn)

        # Progress bar
        self.t1_prog = QProgressBar()
        self.t1_prog.setValue(0)

        # Add everything to tab 1
        self.tab1_layout.addWidget(self.t1_input_btn)
        self.tab1_layout.addLayout(self.t1_sub_1)
        self.tab1_layout.addLayout(self.t1_sub_2)
        self.tab1_layout.addLayout(self.t1_sub_3)
        self.tab1_layout.addLayout(self.t1_sub_4)
        self.tab1_layout.addLayout(self.t1_sub_5)
        self.tab1_layout.addWidget(self.t1_prog)

        # ##### Tab 2 Contents ##### #
        # Input select buttons
        self.t2_sub_1 = QHBoxLayout()
        self.t2_input_btn = QPushButton("Select Input Data")
        self.t2_model_btn = QPushButton("Select Model")
        self.t2_input_btn.clicked.connect(self.sel_t2_file)
        self.t2_model_btn.clicked.connect(self.sel_t2_model)
        self.t2_sub_1.addWidget(self.t2_input_btn)
        self.t2_sub_1.addWidget(self.t2_model_btn)

        # Input and model labels and display
        self.t2_sub_2 = QVBoxLayout()
        self.t2_sub_2_1 = QHBoxLayout()
        self.t2_input_lab = QLabel("Chosen Input:")
        self.t2_input_edit = QLineEdit()
        self.t2_input_edit.setPlaceholderText("Select Input...")
        self.t2_input_edit.setReadOnly(True)
        self.t2_sub_2_1.addWidget(self.t2_input_lab)
        self.t2_sub_2_1.addWidget(self.t2_input_edit)

        self.t2_sub_2_2 = QHBoxLayout()
        self.t2_model_lab = QLabel("Chosen Model:")
        self.t2_model_edit = QLineEdit()
        self.t2_model_edit.setPlaceholderText("Select Model...")
        self.t2_model_edit.setReadOnly(True)
        self.t2_sub_2_2.addWidget(self.t2_model_lab)
        self.t2_sub_2_2.addWidget(self.t2_model_edit)

        self.t2_sub_2.addLayout(self.t2_sub_2_1)
        self.t2_sub_2.addLayout(self.t2_sub_2_2)

        # Metadata and statistics display
        self.t2_sub_3 = QHBoxLayout()
        self.t2_sub_3_1 = QVBoxLayout()

        self.t2_size_lab = QLabel("Blocksize:")
        self.t2_slide_lab = QLabel("Sliding window?:")
        self.t2_err_lab = QLabel("Error rate:")
        self.t2_sub_3_1.addWidget(self.t2_size_lab)
        self.t2_sub_3_1.addWidget(self.t2_slide_lab)
        self.t2_sub_3_1.addWidget(self.t2_err_lab)

        self.t2_stats_edit = QTextEdit()
        self.t2_stats_edit.setText(set_t2_stats_edit())

        # Create widget for label layout, so we can set minimum width (for improved visibility)
        self.t2_sub_3_1_widget = QWidget()
        self.t2_sub_3_1_widget.setLayout(self.t2_sub_3_1)
        self.t2_sub_3_1_widget.setMinimumWidth(int(self.tab2.width() / 2))
        self.t2_sub_3.addWidget(self.t2_sub_3_1_widget)
        self.t2_sub_3.addWidget(self.t2_stats_edit)

        # Functionality buttons
        self.t2_sub_4 = QHBoxLayout()
        self.t2_rst_btn = QPushButton("Reset")
        self.t2_run_btn = QPushButton("Run")
        self.t2_wrt_btn = QPushButton("Write")
        self.t2_gen_btn = QPushButton("Generate")

        self.t2_rst_btn.clicked.connect(self.t2_reset)
        self.t2_run_btn.clicked.connect(self.analyse)
        self.t2_wrt_btn.clicked.connect(self.write_results)

        self.t2_gen_btn.setEnabled(False)

        self.t2_sub_4.addWidget(self.t2_rst_btn)
        self.t2_sub_4.addWidget(self.t2_run_btn)
        self.t2_sub_4.addWidget(self.t2_wrt_btn)
        self.t2_sub_4.addWidget(self.t2_gen_btn)

        # Add everything to tab 2
        self.tab2_layout.addLayout(self.t2_sub_1)
        self.tab2_layout.addLayout(self.t2_sub_2)
        self.tab2_layout.addLayout(self.t2_sub_3)
        self.tab2_layout.addLayout(self.t2_sub_4)

        # ##### Add tabs to main layout ##### #
        self.tabwidget.addTab(self.tab1, "Calculator")
        self.tabwidget.addTab(self.tab2, "Analyser")
        self.main_layout.addWidget(self.tabwidget)

        # ##### Internal variables ##### #
        # Strings representing selected files
        self.t1_file = ""
        self.t2_file = ""
        self.t2_model = ""

        # Chi-square distributions
        self.exp = []
        self.obs = []

    def get_exp(self):
        return self.exp

    def get_obs(self):
        return self.obs

    def set_exp(self, new_exp):
        self.exp = new_exp

    def set_obs(self, new_obs):
        self.obs = new_obs

    def get_t1_file(self):
        return self.t1_file

    def get_t2_file(self):
        return self.t2_file

    def get_t2_model(self):
        return self.t2_model

    def get_t1_slide(self):
        return self.t1_slide_chk.isChecked()

    def get_t1_err(self):
        return self.t1_err_edit.text()

    def get_t1_size(self):
        return self.t1_size_edit.text()

    def set_t1_file(self, name):
        self.t1_file = name
        self.t1_input_edit.setText(name)

    def set_t2_file(self, name):
        self.t2_file = name
        self.t2_input_edit.setText(name)

    def set_t2_model(self, name):
        self.t2_model = name
        self.t2_model_edit.setText(name)

    def sel_t1_file(self):
        file = self.get_file()
        self.set_t1_file(file)

    def sel_t2_file(self):
        file = self.get_file()
        self.set_t2_file(file)

    def sel_t2_model(self):
        model = self.get_file()
        self.set_t2_model(model)

    def set_t2_size(self, value=""):
        self.t2_size_lab.setText("Blocksize: %s" % value)

    def set_t2_slide(self, value=""):
        self.t2_slide_lab.setText("Sliding window?: %s" % value)

    def set_t2_err(self, value=""):
        self.t2_err_lab.setText("Error rate: %s" % value)

    def t1_reset(self):
        """
        Reset tab 1 labels to blank
        :return: None
        """
        self.set_t1_file("")
        self.t1_size_edit.setText("")
        self.t1_err_edit.setText("")

    def t2_reset(self):
        """
        Reset tab 2 labels and results to blank
        :return: None
        """
        self.set_t2_file("")
        self.set_t2_model("")
        self.set_t2_size()
        self.set_t2_slide()
        self.set_t2_err()
        self.t2_stats_edit.setText(set_t2_stats_edit())

    def analyse(self):
        """
        Perform statistical analysis over the chosen BitReps measurements file
        :return: None
        """
        # Obtain metadata
        meta_data = get_meta_data(self.get_t2_file())
        bs, sw, er, nb, oh, avger = meta_data[0], meta_data[1], meta_data[2], meta_data[3], meta_data[4], meta_data[5]

        # Set metadata labels
        self.set_t2_size(bs)
        self.set_t2_slide(sw)
        self.set_t2_err(er)

        # Obtain expected false positives, duplicates, highest rep, ratio
        exp_fps = get_exp_fps(nb, avger)
        exp_dupes = get_exp_dupes(nb, bs)
        highest_rep = get_highest_rep(self.get_t2_file())
        ratio = get_ratio(oh, (exp_fps + exp_dupes))

        # Obtain and unpack chi-square-related information
        chi = calc_chi(self.get_t2_file(), self.get_t2_model())
        chi_val = chi[0]
        self.set_obs(chi[1])
        self.set_exp(chi[2])

        # Write analysis output to display
        self.t2_stats_edit.setText(set_t2_stats_edit(
            str(self.get_exp()),
            str(self.get_obs()),
            round(chi_val, 2),
            exp_fps,
            exp_dupes,
            oh,
            ratio,
            highest_rep
        ))

    def write_results(self):
        """
        Write the results of automated analysis to a .txt file
        :return: None
        """
        results_path = os.path.join(RESULTS_DIR, Path(self.get_t2_file()).stem)
        with open("%s.txt" % results_path, "w") as f:
            f.write(self.t2_stats_edit.toPlainText())

    def measure(self):
        """
        Perform BitReps measurements for the selected file using the user-specified parameters
        :return: None
        """
        bitreps_measure(
            self.get_t1_file(),
            int(self.get_t1_size()),
            self.t1_prog,
            self.get_t1_slide(),
            float(self.get_t1_err())
        )

    def get_file(self):
        """
        Open a file explorer in the current working directory
        :return: Path of selected file
        """
        response = QFileDialog.getOpenFileName(
            parent=self,
            caption="Select input",
            directory=os.getcwd()
        )
        return response[0]


if __name__ == "__main__":
    dir_setup()
    app = QApplication(sys.argv)
    bitReps = BitReps()
    bitReps.show()
    sys.exit(app.exec_())
