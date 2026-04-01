# ==============================================================================
# Functional Data Analysis — EEG Project Makefile
#
# Usage:
#   make all           Run the full pipeline from scratch
#   make deps          Install Python and R dependencies
#   make data          Acquire raw EEG data via datalad
#   make stimuli       Raw EEG → per-stimulus CSVs (~40 min)
#   make assemble      CSVs → subject matrix CSV (~30 sec)
#   make eda           Smoothing + full EDA (~2 min)
#   make presentation  Compile LaTeX slides
#   make clean         Remove EDA outputs and presentation build files
#   make distclean     clean + remove all generated data folders
#   make help          Show all available targets
#
# Prerequisites:
#   Linux/macOS: make is pre-installed
#   Windows:     install via Rtools4 or: choco install make
#                run from Git Bash or Rtools terminal, not PowerShell
#
# ==============================================================================

# --- Configuration -----------------------------------------------------------

PYTHON     := python3
RSCRIPT    := Rscript
LATEXMK    := latexmk

NOTEBOOKS  := Notebooks
EDA_DIR    := EDA
PRES_DIR   := Presentations
OUT_DIR    := $(EDA_DIR)/outputs

# Input/output files
SUBJECT_CSV   := $(EDA_DIR)/Flanker_stimulus_FC1_channel.csv
SUBJECT_META  := $(EDA_DIR)/subject_metadata.csv
FD_SMOOTH     := $(OUT_DIR)/fd_smooth.rds
STIMULI_DIR   := ds006018_per_stimuli
FUNC_DIR      := ds006018_functional
RAW_DATA_DIR  := ds006018
PRESENTATION  := $(PRES_DIR)/main.pdf

# --- Phony targets -----------------------------------------------------------

.PHONY: all deps deps-python deps-r data stimuli assemble eda presentation \
        clean distclean help

# --- Default: full pipeline --------------------------------------------------

all: deps data stimuli assemble eda presentation
	@echo ""
	@echo "=== Full pipeline complete. ==="

# --- Help --------------------------------------------------------------------

help:
	@echo ""
	@echo "  make all           Run the full pipeline from scratch"
	@echo "  make deps          Install Python and R dependencies"
	@echo "  make data          Acquire raw EEG data via datalad"
	@echo "  make stimuli       Raw EEG -> per-stimulus CSVs (~40 min)"
	@echo "  make assemble      CSVs -> subject matrix CSV (~30 sec)"
	@echo "  make eda           Smoothing + full EDA (~2 min)"
	@echo "  make presentation  Compile LaTeX slides"
	@echo "  make clean         Remove EDA outputs and presentation build files"
	@echo "  make distclean     clean + remove all generated data folders"
	@echo ""

# --- Step 1: Dependencies ----------------------------------------------------

deps: deps-python deps-r

deps-python:
	$(PYTHON) -m venv .venv
	.venv/bin/pip install -r requirements.txt || .venv/Scripts/pip install -r requirements.txt
	@echo "Python dependencies installed."

deps-r:
	$(RSCRIPT) -e "if (!requireNamespace('renv', quietly=TRUE)) install.packages('renv'); renv::restore()"
	@echo "R dependencies installed."

# --- Step 2: Raw data --------------------------------------------------------

data: $(RAW_DATA_DIR)/.datalad

$(RAW_DATA_DIR)/.datalad:
	git submodule update --init
	cd $(RAW_DATA_DIR) && datalad get .
	@echo "Raw EEG data acquired."

# --- Step 3: Per-stimulus CSVs (Python) --------------------------------------

stimuli: $(STIMULI_DIR)/.done

$(STIMULI_DIR)/.done: $(RAW_DATA_DIR)/.datalad
	cd $(NOTEBOOKS) && $(PYTHON) -m jupyter nbconvert --to notebook --execute 03_data_preparation.ipynb --output /dev/null
	@touch $@
	@echo "Per-stimulus CSVs generated."

# --- Step 4: Assemble subject matrix (R) -------------------------------------

assemble: $(SUBJECT_CSV)

$(SUBJECT_CSV): $(EDA_DIR)/assemble_subject_flanker_S2_FC1.R $(STIMULI_DIR)/.done
	cd $(EDA_DIR) && $(RSCRIPT) assemble_subject_flanker_S2_FC1.R
	@echo "Subject matrix assembled: $(SUBJECT_CSV)"

# --- Step 5: Smoothing + EDA (R) ---------------------------------------------

eda: $(FD_SMOOTH)

$(FD_SMOOTH): $(EDA_DIR)/Smoothing_and_EDA.R $(SUBJECT_CSV)
	cd $(EDA_DIR) && $(RSCRIPT) Smoothing_and_EDA.R
	@echo "EDA complete. Outputs in $(OUT_DIR)/"

# --- Step 6: LaTeX presentation ----------------------------------------------

presentation: $(PRESENTATION)

$(PRESENTATION): $(PRES_DIR)/main.tex $(FD_SMOOTH)
	cp $(OUT_DIR)/*.pdf $(PRES_DIR)/ 2>/dev/null || true
	cd $(PRES_DIR) && $(LATEXMK) -xelatex -interaction=nonstopmode main.tex
	@echo "Presentation compiled: $(PRESENTATION)"

# --- Clean -------------------------------------------------------------------

clean:
	rm -f $(OUT_DIR)/*.pdf $(OUT_DIR)/*.rds $(OUT_DIR)/*.txt
	rm -f $(SUBJECT_CSV) $(SUBJECT_META)
	rm -f $(PRES_DIR)/*.aux $(PRES_DIR)/*.log $(PRES_DIR)/*.nav
	rm -f $(PRES_DIR)/*.out $(PRES_DIR)/*.snm $(PRES_DIR)/*.toc
	rm -f $(PRES_DIR)/*.fls $(PRES_DIR)/*.fdb_latexmk $(PRES_DIR)/*.xdv
	rm -f $(PRES_DIR)/main.pdf
	@echo "Cleaned EDA outputs and presentation build files."

distclean: clean
	rm -rf $(STIMULI_DIR)
	rm -rf $(FUNC_DIR)
	@echo "Also removed per-stimulus CSVs and functional .rds files."
