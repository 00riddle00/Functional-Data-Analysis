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

NOTEBOOKS   := Notebooks
EDA_DIR     := EDA
HT_DIR      := HT
PRES_DIR    := Presentations
EDA_OUT_DIR := $(EDA_DIR)/outputs
HT_OUT_DIR  := $(HT_DIR)/outputs

RAW_DATA_URL := https://github.com/OpenNeuroDatasets/ds006018.git

# Input/output files
SUBJECT_CSV     := $(EDA_DIR)/Flanker_stimulus_FC1_channel.csv
SUBJECT_META    := $(EDA_DIR)/subject_metadata.csv
FD_SMOOTH       := $(EDA_OUT_DIR)/fd_smooth.rds
# HT_01 also tracks HT_02, HT_03, and HT_04 since they are generated together in the same script.
HT_PLOTS        := $(HT_OUT_DIR)/HT_01_group_comparison.pdf
# HT_05 also tracks HT_06, HT_07, and HT_08.
HT_PLOTS_GENDER := $(HT_OUT_DIR)/HT_05_gender_group_comparison.pdf
# HT_09 also tracks HT_10, HT_11, and HT_12.
HT_PLOTS_SES    := $(HT_OUT_DIR)/HT_09_ses_group_comparison.pdf
STIMULI_DIR     := ds006018_per_stimuli
FUNC_DIR        := ds006018_functional
RAW_DATA_DIR    := ds006018
PRESENTATION_1  := $(PRES_DIR)/presentation_1st.pdf
PRESENTATION_2  := $(PRES_DIR)/presentation_2nd.pdf
REG_DIR         := REG
# REG_01 also tracks REG_02, REG_03, and REG_04.
REG_PLOTS       := $(REG_DIR)/outputs/REG_01_coefficients.pdf
REPORTS_DIR     := Reports
REPORT          := $(REPORTS_DIR)/final_report.pdf

# --- Phony targets -----------------------------------------------------------

.PHONY: all deps deps-python deps-r sync-requirements sync-uv data stimuli functional assemble \
	eda hypothesis_testing presentation_1 presentation_2 regression report clean distclean help

# --- Default: full pipeline --------------------------------------------------

all: deps data stimuli assemble eda hypothesis_testing presentation_1 presentation_2
	@echo ""
	@echo "=== Full pipeline complete. ==="

# --- Help --------------------------------------------------------------------

help:
	@echo ""
	@echo "  make all                  Run the full pipeline from scratch"
	@echo "  make deps                 Install Python and R dependencies"
	@echo "  make export-requirements  Export uv lockfile to requirements.txt"
	@echo "  make import-requirements  Import requirements.txt changes into uv"
	@echo "  make data                 Acquire raw EEG data via datalad"
	@echo "  make stimuli              Raw EEG -> per-stimulus CSVs (~40 min)"
	@echo "  make functional           Generate F7 .rds files (optional, ~40 min)"
	@echo "  make assemble             CSVs -> subject matrix CSV"
	@echo "  make eda                  Smoothing + full EDA"
	@echo "  make hypothesis_testing   Run all hypothesis testing scripts"
	@echo "  make presentation_1       Compile LaTeX slides for 1st presentation"
	@echo "  make presentation_2       Compile LaTeX slides for 2nd presentation"
	@echo "  make regression           Run regression analyses"
	@echo "  make report               Compile final LaTeX report PDF"
	@echo "  make clean                Remove EDA outputs and presentation build files"
	@echo "  make distclean            Clean + remove all generated data folders (caution)"
	@echo "  make clean-env            Remove Python venv and R library (for testing)"
	@echo "  make help                 Show this message"
	@echo ""

# --- Step 1: Dependencies ----------------------------------------------------

deps: deps-python deps-r

deps-python: .venv/.stamp

deps-r: renv/.stamp

# Rebuilds venv if requirements.txt has a newer timestamp than .stamp.
# Note: git operations (pull, checkout) can update file timestamps,
# causing an unnecessary rebuild. This is harmless but slow —
# `python -m venv .venv` reinitializes the venv structure without wiping
# installed packages, and pip skips already-installed dependencies.
.venv/.stamp: requirements.txt
	$(PYTHON) -m venv .venv
	$(VENV_BIN)/pip install --upgrade pip
	$(VENV_BIN)/pip install -r requirements.txt
	$(STAMP_DATE) > $@
	@echo "Python dependencies installed."

# Rebuilds R library if renv.lock has a newer timestamp than .stamp.
# Note: git operations (pull, checkout) can update file timestamps,
# causing an unnecessary rebuild. This is harmless but slow —
# renv::restore() checks each package and skips already-installed ones.
# The renv/ directory itself is not wiped; only missing or outdated
# packages are reinstalled.
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

# Export uv-managed dependencies to classic requirements.txt
export-requirements:
	uv export \
	  --format requirements-txt \
	  --no-hashes \
	  --no-header \
	  --no-annotate \
	  > requirements.txt

# Import updated requirements.txt into uv workflow
import-requirements:
	uv add -r requirements.txt
	uv lock
	uv sync

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
#
# NOTE: The notebook file (03_data_preparation.ipynb) is intentionally NOT
# listed as a dependency. Editing the notebook will not trigger a re-run.
# This avoids accidental 40+ minute re-runs from minor notebook edits
# (e.g., adding comments or formatting). To force a re-run after a
# meaningful notebook change, delete the stamp file:
#   rm ds006018_per_stimuli/.stamp && make stimuli
stimuli: $(STIMULI_DIR)/.stamp

# TODO: consider replacing notebook execution with .py/.R scripts for pipeline
$(STIMULI_DIR)/.stamp: $(RAW_DATA_DIR)/.stamp
	$(VENV_BIN)/jupyter nbconvert \
		--to notebook \
		--execute \
		--inplace \
		$(NOTEBOOKS)/03_data_preparation.ipynb
	@mkdir -p $(dir $@)
	$(STAMP_DATE) > $@
	@echo "Per-stimulus CSVs generated."

# --- Step 3b (optional): F7 .rds files ---------------------------------------
#
# NOTE: Same as above - notebook edits do not trigger re-runs.
# To force: rm ds006018_functional/.stamp && make functional
functional: $(FUNC_DIR)/.stamp

# TODO: consider replacing notebook execution with .py/.R scripts for pipeline
$(FUNC_DIR)/.stamp: $(STIMULI_DIR)/.stamp
	$(VENV_BIN)/jupyter nbconvert \
		--to notebook \
		--execute \
		--inplace \
		$(NOTEBOOKS)/04_data_preparation_R.ipynb
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
	@echo "EDA complete. Outputs in $(EDA_OUT_DIR)/"

# --- Step 6: Hypothesis testing ----------------------------------------------

hypothesis_testing: $(HT_PLOTS) $(HT_PLOTS_GENDER) $(HT_PLOTS_SES)

$(HT_PLOTS): $(HT_DIR)/Hypothesis_testing_ADHD.R $(FD_SMOOTH)
	$(RSCRIPT) $(HT_DIR)/Hypothesis_testing_ADHD.R
	@echo "Hypothesis testing (ADHD) complete. Outputs in $(HT_OUT_DIR)/"

$(HT_PLOTS_GENDER): $(HT_DIR)/Hypothesis_testing_Gender.R $(FD_SMOOTH)
	$(RSCRIPT) $(HT_DIR)/Hypothesis_testing_Gender.R
	@echo "Hypothesis testing (Gender) complete. Outputs in $(HT_OUT_DIR)/"

$(HT_PLOTS_SES): $(HT_DIR)/Hypothesis_testing_SES.R $(FD_SMOOTH)
	$(RSCRIPT) $(HT_DIR)/Hypothesis_testing_SES.R
	@echo "Hypothesis testing (SES) complete. Outputs in $(HT_OUT_DIR)/"

# --- Step 7: LaTeX presentations ---------------------------------------------

presentation_1: $(PRESENTATION_1)

$(PRESENTATION_1): $(PRES_DIR)/presentation_1st.tex $(FD_SMOOTH)
	cp $(EDA_OUT_DIR)/*.pdf $(PRES_DIR)/ 2>/dev/null || true
	$(LATEXMK) -xelatex -interaction=nonstopmode -outdir=$(PRES_DIR) $(PRES_DIR)/presentation_1st.tex
	@echo "Presentation compiled: $(PRESENTATION_1)"

presentation_2: $(PRESENTATION_2)

$(PRESENTATION_2): $(PRES_DIR)/presentation_2nd.tex $(HT_PLOTS)
	cp $(HT_OUT_DIR)/*.pdf $(PRES_DIR)/ 2>/dev/null || true
	$(LATEXMK) -xelatex -interaction=nonstopmode -outdir=$(PRES_DIR) $(PRES_DIR)/presentation_2nd.tex
	@echo "Presentation compiled: $(PRESENTATION_2)"

# --- Step 8: Regression ------------------------------------------------------

regression: $(REG_PLOTS)

$(REG_PLOTS): $(REG_DIR)/FDA_regression.R $(FD_SMOOTH)
	$(RSCRIPT) $(REG_DIR)/FDA_regression.R
	@echo "Regression complete. Outputs in $(REG_DIR)/outputs/"

# --- Step 9: Final report ----------------------------------------------------

report: $(REPORT)

$(REPORT): $(REPORTS_DIR)/final_report.tex $(HT_PLOTS) $(HT_PLOTS_GENDER) $(HT_PLOTS_SES) $(REG_PLOTS)
	cp $(EDA_OUT_DIR)/*.pdf $(REPORTS_DIR)/ 2>/dev/null || true
	cp $(HT_OUT_DIR)/*.pdf $(REPORTS_DIR)/ 2>/dev/null || true
	cp $(REG_DIR)/outputs/*.pdf $(REPORTS_DIR)/ 2>/dev/null || true
	$(LATEXMK) -pdf -interaction=nonstopmode -outdir=$(REPORTS_DIR) $(REPORTS_DIR)/final_report.tex
	@echo "Report compiled: $(REPORT)"

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
	rm -f $(EDA_OUT_DIR)/*.pdf $(EDA_OUT_DIR)/*.rds $(EDA_OUT_DIR)/*.txt
	rm -f $(HT_OUT_DIR)/*.pdf $(HT_OUT_DIR)/*.rds $(HT_OUT_DIR)/*.txt
	rm -f $(REG_DIR)/outputs/*.pdf
	rm -f $(REPORTS_DIR)/*.pdf
	rm -f $(SUBJECT_CSV) $(SUBJECT_META)
	git clean -fdX -- $(PRES_DIR)
	@echo "Cleaned EDA, HT outputs and presentation build files."

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
