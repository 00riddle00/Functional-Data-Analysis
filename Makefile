# vim: set ft=make tw=100 nu noet ts=8 sw=8:
# =============================================================================
# Functional Data Analysis — EEG Project Makefile
#
# Prerequisites:
#   Linux/macOS: make is pre-installed
#   Windows:     install via Rtools4 or: choco install make
#                run from Git Bash or Rtools terminal, not PowerShell
# =============================================================================

# --- Configuration -----------------------------------------------------------

PYTHON     := python3
RSCRIPT    := Rscript
LATEXMK    := latexmk
STAMP_DATE := date "+%F %T %Z"

ifeq ($(OS),Windows_NT)
	VENV_BIN := $(CURDIR)/.venv/Scripts
else
	VENV_BIN := $(CURDIR)/.venv/bin
endif

NOTEBOOKS  := Notebooks
EDA_DIR    := EDA
PRES_DIR   := Presentations
OUT_DIR    := $(EDA_DIR)/outputs

RAW_DATA_URL := https://github.com/OpenNeuroDatasets/ds006018.git

# Input/output files
SUBJECT_CSV   := $(EDA_DIR)/Flanker_stimulus_FC1_channel.csv
SUBJECT_META  := $(EDA_DIR)/subject_metadata.csv
FD_SMOOTH     := $(OUT_DIR)/fd_smooth.rds
STIMULI_DIR   := ds006018_per_stimuli
FUNC_DIR      := ds006018_functional
RAW_DATA_DIR  := ds006018
PRESENTATION  := $(PRES_DIR)/presentation_1st.pdf

# --- Phony targets -----------------------------------------------------------

.PHONY: all deps deps-python deps-r data stimuli functional assemble eda presentation clean \
	distclean help

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
	@echo "  make functional    Generate F7 .rds files (optional, ~40 min)"
	@echo "  make assemble      CSVs -> subject matrix CSV"
	@echo "  make eda           Smoothing + full EDA"
	@echo "  make presentation  Compile LaTeX slides"
	@echo "  make clean         Remove EDA outputs and presentation build files"
	@echo "  make distclean     clean + remove all generated data folders (caution)"
	@echo "  make clean-env     Remove Python venv and R library (for testing)"
	@echo "  make help          show this message"
	@echo ""

# --- Step 1: Dependencies ----------------------------------------------------

deps: deps-python deps-r

deps-python: .venv/.stamp

deps-r: renv/.stamp

.venv/.stamp: requirements.txt
	$(PYTHON) -m venv .venv
	$(VENV_BIN)/pip install --upgrade pip
	$(VENV_BIN)/pip install -r requirements.txt
	$(STAMP_DATE) > $@
	@echo "Python dependencies installed."

renv/.stamp: renv.lock
	@mkdir -p renv
	$(RSCRIPT) -e "\
	  if (!requireNamespace('renv', quietly=TRUE)) install.packages('renv'); \
	  renv::restore(prompt = FALSE); \
	  if (!requireNamespace('IRkernel', quietly=TRUE)) install.packages('IRkernel'); \
	  IRkernel::installspec( \
	    name = paste0('r-', basename(getwd())), \
	    displayname = paste0('R (', basename(getwd()), ')') \
	  ); \
	"
	$(STAMP_DATE) > $@
	@echo "R dependencies installed and IRkernel registered."

# --- Step 2: Raw data --------------------------------------------------------
#
# .stamp is a stamp/sentinel file (= empty target) created after datalad
# downloads the full dataset (~10 GB). Without it, make cannot distinguish
# between "repo cloned but data not fetched" and "data fully fetched," since
# .datalad/ exists in both cases. It also contains a timestamp of the last
# successful data acquisition, which can be useful for debugging and tracking
# when the data was last updated.
data: $(RAW_DATA_DIR)/.stamp

$(RAW_DATA_DIR)/.stamp:
	@if [ ! -d "$(RAW_DATA_DIR)" ]; then \
		git clone $(RAW_DATA_URL) $(RAW_DATA_DIR); \
	fi
	$(VENV_BIN)/datalad get -d $(RAW_DATA_DIR) $(RAW_DATA_DIR)
	$(STAMP_DATE) > $@
	@echo "Raw EEG data acquired."

# --- Step 3: Per-stimulus CSVs -----------------------------------------------

stimuli: $(STIMULI_DIR)/.stamp

$(STIMULI_DIR)/.stamp: $(RAW_DATA_DIR)/.stamp
	cd $(NOTEBOOKS) && $(VENV_BIN)/jupyter nbconvert \
		--to notebook \
		--execute \
		--inplace \
		03_data_preparation.ipynb
	$(STAMP_DATE) > $@
	@echo "Per-stimulus CSVs generated."

# TODO: consider replacing notebook execution with .py/.R scripts for pipeline
#$(STIMULI_DIR)/.stamp: $(RAW_DATA_DIR)/.stamp
	#$(VENV_BIN)/jupyter nbconvert \
		#--to notebook \
		#--execute \
		#--inplace \
		#$(NOTEBOOKS)/03_data_preparation.ipynb
	#@mkdir -p $(dir $@)
	#$(STAMP_DATE) > $@
	#@echo "Per-stimulus CSVs generated."

# --- Step 3b (optional): F7 .rds files ---------------------------------------

functional: $(FUNC_DIR)/.stamp

# TODO: consider replacing notebook execution with .py/.R scripts for pipeline
$(FUNC_DIR)/.stamp: $(STIMULI_DIR)/.stamp
	$(VENV_BIN)/jupyter nbconvert \
		--to notebook \
		--execute $(NOTEBOOKS)/04_data_preparation_R.ipynb \
		--output /dev/null
	$(STAMP_DATE) > $@
	@echo "Functional .rds files generated."

# --- Step 4: Assemble subject matrix -----------------------------------------

assemble: $(SUBJECT_CSV)

$(SUBJECT_CSV): $(EDA_DIR)/assemble_subject_flanker_S2_FC1.R $(STIMULI_DIR)/.stamp
	$(RSCRIPT) $(EDA_DIR)/assemble_subject_flanker_S2_FC1.R
	@echo "Subject matrix assembled: $(SUBJECT_CSV)"

# --- Step 5: Smoothing + EDA -------------------------------------------------

eda: $(FD_SMOOTH)

$(FD_SMOOTH): $(EDA_DIR)/Smoothing_and_EDA.R $(SUBJECT_CSV)
	$(RSCRIPT) $(EDA_DIR)/Smoothing_and_EDA.R
	@echo "EDA complete. Outputs in $(OUT_DIR)/"

# --- Step 6: LaTeX presentation ----------------------------------------------

presentation: $(PRESENTATION)

$(PRESENTATION): $(PRES_DIR)/presentation_1st.tex $(FD_SMOOTH)
	cp $(OUT_DIR)/*.pdf $(PRES_DIR)/ 2>/dev/null || true
	$(LATEXMK) -xelatex -interaction=nonstopmode -outdir=$(PRES_DIR) $(PRES_DIR)/presentation_1st.tex
	@echo "Presentation compiled: $(PRESENTATION)"

# --- Clean -------------------------------------------------------------------

# TODO:
# `make clean` currently removes files that are tracked by Git (EDA outputs, CSVs, etc.).
# This is intentional for now, but results in `git status` showing deletions.
# In the future, consider:
#   - moving generated artifacts to .gitignore, OR
#   - keeping only a curated subset of outputs in Git, OR
#   - using a tool like DVC for data/artifact versioning.
# Current workflow: run `git restore EDA/` after `make clean` if needed.
clean:
	rm -f $(OUT_DIR)/*.pdf $(OUT_DIR)/*.rds $(OUT_DIR)/*.txt
	rm -f $(SUBJECT_CSV) $(SUBJECT_META)
	git clean -fdX -- $(PRES_DIR)
	@echo "Cleaned EDA outputs and presentation build files."

# Caution: will remove all generated data folders, which are not tracked by Git.
distclean: clean
	rm -rf $(STIMULI_DIR)
	rm -rf $(FUNC_DIR)
	rm -f $(RAW_DATA_DIR)/.stamp
	@echo "Also removed generated data folders and sentinel files."

# Completely remove Python venv and R library (for testing purposes).
clean-env:
	rm -rf .venv
	rm -rf renv/library
	rm -f renv/.stamp
	@echo "Removed Python venv, R library, and dependency stamps."
