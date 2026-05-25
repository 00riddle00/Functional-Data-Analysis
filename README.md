<!-- vim: set ft=markdown tw=88 nu ai et ts=2 sw=2: -->

# Functional Data Analysis of EEG Data

**Course:** Functional Data Analysis (10 ECTS), VU MIF, Spring 2026

**Study program:** MSc Data Science

**Team:** [Jonas Adomaitis](https://github.com/JonasIBM), [Gedas
Beržinskas](https://github.com/Berzinskass), [Tomas
Giedraitis](https://github.com/00riddle00)

Applying Functional Data Analysis methods to EEG brain recordings from 127 young adults
performing a Flanker cognitive task. Data source: [OpenNeuro
ds006018](https://github.com/OpenNeuroDatasets/ds006018) — *Cognitive Electrophysiology
in Socioeconomic Context in Adulthood* (Isbell et al., 2025,
[doi:10.1038/s41597-025-05209-z](https://www.nature.com/articles/s41597-025-05209-z)).

- [Quick start](#quick-start)
- [Prerequisites](#prerequisites)
- [Setup](#setup)
- [Pipeline](#pipeline)
- [Folder structure](#folder-structure)
- [R / Python Pipeline](#r--python-pipeline)
- [Notes](#notes)
- [Make targets](#make-targets)

## Quick start

```bash
$ git clone https://github.com/00riddle00/Functional-Data-Analysis
$ cd Functional-Data-Analysis
$ git clone https://github.com/OpenNeuroDatasets/ds006018
$ make deps
$ make all
```

For individual steps, run `$ make help`.

---

## Prerequisites

### Git

- **Linux:** `$ sudo apt install git` (Ubuntu/Debian) or use your distro's package
  manager
- **macOS:** `$ brew install git`
- **Windows:** https://git-scm.com/download/win — after install, run: `$ git config
  --global core.longpaths true`

### R (4.5+)

- **Linux:** https://cran.r-project.org/bin/linux/ (follow distro-specific instructions)
- **macOS:** https://cran.r-project.org/bin/macosx/
- **Windows:** https://cran.r-project.org/bin/windows/base/

### Python (3.14+)

- **Linux:** `$ sudo apt install python3 python3-venv python3-pip` (or use your distro's
  package manager)
- **macOS:** `$ brew install python`
- **Windows:** https://www.python.org — check "Add to PATH" during install

### DataLad and git-annex

Needed to acquire the raw 10 GB EEG data.

- **Linux:** `$ sudo apt install datalad git-annex`
- **macOS:** `(.venv) $ brew install git-annex && pip install datalad`
- **Windows:** Install git-annex from https://git-annex.branchable.com/install/Windows/
  then `(.venv) $ pip install datalad`.

### LaTeX

Needed to compile the presentations and final report locally.

- **Linux:** `$ sudo apt install texlive-full` (or minimal:
  `$ sudo apt install texlive-base texlive-latex-recommended texlive-latex-extra
  texlive-fonts-recommended texlive-xetex latexmk`)
- **macOS:** https://www.tug.org/mactex/ or `$ brew install --cask mactex`
- **Windows:** https://miktex.org/download (select "Install missing packages on the
  fly"). Alternatively, TeX Live: https://www.tug.org/texlive/windows.html

### System libraries (Linux only)

Some R packages require system libraries. Install before running `> renv::restore()` in
R console:

```bash
$ sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
  libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```

### Jupyter R kernel

Needed to run R-based notebooks in Jupyter. Run in R console:

```r
> install.packages("IRkernel")
IRkernel::installspec(
  name = paste0("r-", basename(getwd())),
  displayname = paste0("R (", basename(getwd()), ")")
)
```

This creates a project-specific kernel name. Renv activation is handled by the first
cell in each R notebook.

### Make

Needed to use the Makefile for automation.

- **Linux:** Usually pre-installed. If not: `$ sudo apt install make`
- **macOS:** Included with Xcode command line tools (`$ xcode-select --install`) or just
  `$ brew install make`
- **Windows:** Comes with Rtools4 (which you likely have for R package compilation).
  Otherwise: `$ choco install make`. Run from Git Bash or Rtools terminal, not
  PowerShell.

---

## Setup

### Clone

```bash
$ git clone https://github.com/00riddle00/Functional-Data-Analysis
$ cd Functional-Data-Analysis
$ git clone https://github.com/OpenNeuroDatasets/ds006018
```

### Install dependencies

**R:**

Run in R console:

```r
> install.packages("renv")
> renv::restore()
```

**Python:**

```bash
$ python3 -m venv .venv
$ source .venv/bin/activate      # Linux/macOS
# $ .\venv\Scripts\activate      # Windows
(.venv) $ pip install -r requirements.txt
```

**uv vs. pip:**

This project uses `uv` internally for Python dependency locking. `requirements.txt`
is kept as a `pip`-compatible export for university setup.

After changing dependencies with `uv`, run:

```bash
$ make export-requirements
```

After changing dependencies with pip and updating `requirements.txt`, run:

```bash
$ make import-requirements
```

---

## Important Notes

- `ds006018/` is a cloned dataset repository
  (https://github.com/OpenNeuroDatasets/ds006018.git). It is not tracked by this
  repository and should not be committed.
- `ds006018_per_stimuli/` is in `.gitignore` (54 GB). Regenerate from step 2.
- `ds006018_functional/` is in `.gitignore`. Contains F7-channel .rds files. Currently
  unused — the main analysis uses FC1 extracted directly from the CSVs.

---

### Makefile and stamp files

The Makefile uses sentinel `.stamp` files to track whether long-running steps have
already completed. If you have generated outputs from a previous run but the `.stamp`
files are missing (e.g. after a fresh clone or a `git checkout`), `make` will re-run the
corresponding steps unnecessarily. To prevent this, touch the relevant stamps manually:

```bash
touch ds006018_per_stimuli/.stamp
touch EDA/Flanker_stimulus_FC1_channel.csv
touch EDA/outputs/fd_smooth.rds
touch HT/outputs/HT_01_group_comparison.pdf
touch HT/outputs/HT_05_gender_group_comparison.pdf
touch HT/outputs/HT_09_ses_group_comparison.pdf
touch REG/outputs/REG_01_coefficients.pdf
```

This tells `make` that all pipeline outputs are up to date, so subsequent targets such
as `make presentation_3` or `make report` will skip straight to compilation.

---

## Pipeline

The pipeline has 9 steps. Each depends on the output of the previous one. Run `$ make
all` to execute the full pipeline, or run each step individually.

### Step 1: Acquire raw EEG data (10+ GB, takes a while)

**DataLad:**

```bash
(.venv) $ datalad get -d ds006018 ds006018
```

### Step 2: Generate per-stimulus CSVs (~40 min)

```bash
(.venv) $ jupyter nbconvert \
            --to notebook \
            --execute 03_data_preparation.ipynb \
            --output /dev/null
```

Reads raw BrainVision EEG files, filters, epochs, and exports one CSV per subject per
task per stimulus.

- **Input:** `ds006018/sub-XXX/eeg/*.vhdr`
- **Output:** `ds006018_per_stimuli/sub-XXX/task-XXX_Stimulus_SXX.csv` (50+ GB)

### Step 2b (optional): Generate smoothed .rds files (~40 min)

Only needed if you want F7-channel functional data objects. The main analysis (Steps
3–9) does not use these — it extracts FC1 directly from the CSVs.

```bash
(.venv) $ jupyter nbconvert \
            --to notebook \
            --execute 04_data_preparation_R.ipynb \
            --output /dev/null
```

- **Input:** `ds006018_per_stimuli/sub-XXX/*.csv`
- **Output:** `ds006018_functional/sub-XXX/*.rds` (45+ MB, F7 channel only)

**Note:** Some output filenames contain `:` characters (e.g., `LostSamples:264.rds`)
which are invalid on Windows. The last cell of the notebook renames these automatically.

### Step 3: Assemble the across-subjects data matrix

```bash
$ Rscript EDA/assemble_subject_flanker_S2_FC1.R
```

Extracts FC1 channel from flanker S2 stimulus, averages epochs if needed, combines all
62 subjects into one CSV.

- **Input:** `ds006018_per_stimuli/sub-XXX/task-flanker_Stimulus_S2.csv`
- **Output:** `EDA/Flanker_stimulus_FC1_channel.csv` + `EDA/subject_metadata.csv`

### Step 4: Smoothing + Exploratory Data Analysis (EDA)

```bash
$ Rscript EDA/Smoothing_and_EDA.R
```

Applies B-spline smoothing with GCV-selected λ, then runs the full functional EDA:
mean/SD, covariance, FPCA, depth, outlier detection, boxplots, rainbow plots.

- **Input:** `EDA/Flanker_stimulus_FC1_channel.csv`
- **Output:** 30+ PDF plots + 1 text file + `fd_smooth.rds` in `EDA/outputs/`

### Step 5: Hypothesis testing

There were three hypotheses made:

- **H1:** The shape and amplitude of the functional EEG curves over time are
  significantly associated with **ADHD symptom presence (ADHD vs. non-ADHD
  classification)**.
- **H2:** There is a significant difference in the mean functional trajectories of the
  EEG time-series curves **between genders** across the trial time window.
- **H3:** Self-reported **socioeconomic status (SES)** significantly predicts the
  temporal dynamics and amplitude of the mean functional ERP curves.

```bash
$ Rscript HT/Hypothesis_testing_ADHD.R
$ Rscript HT/Hypothesis_testing_Gender.R
$ Rscript HT/Hypothesis_testing_SES.R
```

Tests whether functional EEG curves differ significantly between groups using pointwise
Z-test, L2-norm test, F-type test, and permutation test (ADHD and gender), and
functional one-way ANOVA (SES).

- **Input:** `EDA/outputs/fd_smooth.rds`
- **Output:** 10 PDF plots in `HT/outputs/`

### Step 6: Regression

```bash
$ Rscript REG/FDA_regression.R
```

Fits a function-on-scalar regression model with ADHD status and gender as predictors.
Also fits a Bayesian alternative using `bayes_fosr`.

- **Input:** `EDA/outputs/fd_smooth.rds`
- **Output:** 4 PDF plots in `REG/outputs/`

### Step 7: Compile presentations

```bash
$ latexmk -xelatex -interaction=nonstopmode -outdir=Presentations Presentations/presentation_1st.tex
$ latexmk -xelatex -interaction=nonstopmode -outdir=Presentations Presentations/presentation_2nd.tex
```

The presentations use Beamer with XeLaTeX. If `latexmk` is unavailable, run
`$ xelatex <file>.tex` twice manually.

### Step 8: Compile final report

```bash
$ latexmk -pdf -interaction=nonstopmode -outdir=Reports Reports/final_report.tex
```

Compiles the JMLR-style final report. Figure assets are copied automatically from
pipeline outputs by `make report`.

- **Input:** `Reports/final_report.tex`, figure PDFs from EDA/HT/REG outputs
- **Output:** `Reports/final_report.pdf`

### LaTeX troubleshooting

- **xeCJK error:** If you don't need CJK (Chinese/Japanese/Korean) support, remove
  `\usepackage{xeCJK}` from the relevant `.tex` file.
- **Font error on Linux:** Install Utopia font: `$ sudo apt install
  texlive-fonts-extra`, or remove `\usepackage{utopia}` from the `.tex` file.
- **Missing packages on Windows (MiKTeX):** MiKTeX installs missing packages
  automatically on first compile. If prompted, click "Install".
- **Missing packages on Linux (TeX Live):** `$ tlmgr install <package-name>`

---

## Analysis Notes

- The analysis uses **channel FC1** (frontal-central, relevant for attentional conflict
  processing in the Flanker task) and **stimulus S2**.
- 62 subjects had flanker S2 data: 55 with 1 trial, 7 with 2 trials (averaged).

---

## Folder structure

```
.
├── ds006018/                                  # Cloned dataset repo — raw EEG data (10+ GB via datalad)
├── ds006018_per_stimuli/                      # .gitignore — intermediate CSVs (50+ GB)
├── ds006018_functional/                       # .gitignore — F7 functional data objects (.rds and .pkl) (unused in main analysis)
├── EDA/
│   ├── assemble_subject_flanker_S2_FC1.R      # Step 3 (R): CSVs → subject matrix
│   ├── assemble_subject_flanker_S2_FC1.ipynb  # Step 3 (Python): CSVs → subject matrix
│   ├── Smoothing_and_EDA.R                    # Step 4 (R): smoothing + full EDA
│   ├── Smoothing_and_EDA.ipynb                # Step 4 (Python): smoothing + full EDA
│   ├── Flanker_stimulus_FC1_channel.csv       # 62 subjects × 501 time points
│   ├── subject_metadata.csv                   # epoch counts per subject
│   └── outputs/                               # 30+ PDFs, fd_smooth.rds, fd_smooth.pkl, text results
├── HT/
│   ├── Hypothesis_testing_ADHD.R              # Step 5a (R): ADHD group comparison
│   ├── Hypothesis_testing_ADHD.ipynb          # Step 5a (Python): ADHD group comparison
│   ├── Hypothesis_testing_Gender.R            # Step 5b (R): gender group comparison
│   ├── Hypothesis_testing_Gender.ipynb        # Step 5b (Python): gender group comparison
│   ├── Hypothesis_testing_SES.R               # Step 5c (R): SES one-way ANOVA
│   ├── Hypothesis_testing_SES.ipynb           # Step 5c (Python): SES one-way ANOVA
│   ├── trace.R                                # helper functions (R)
│   ├── trace.py                               # helper functions (Python)
│   ├── Ztwosample.R                           # helper functions (R)
│   ├── z_two_sample.py                        # helper functions (Python)
│   ├── L2stattwosample.R                      # helper functions (R)
│   ├── l2_stat_two_sample.py                  # helper functions (Python)
│   ├── Fstattwosample.R                       # helper functions (R)
│   ├── f_stat_two_sample.py                   # helper functions (Python)
│   └── outputs/                               # 10 PDF plots
├── REG/
│   ├── FDA_regression.R                       # Step 6 (R): function-on-scalar regression
│   ├── FDA_regression.ipynb                   # Step 6 (Python): function-on-scalar regression
│   └── outputs/                               # 4 PDF plots
├── Reports/
│   ├── final_report.tex                       # LaTeX source
│   ├── final_report.pdf                       # Compiled report
│   ├── bibliography.bib                       # References
│   ├── jmlr2e.sty                             # JMLR style file
│   └── *.pdf                                  # Figure assets copied from pipeline outputs
├── Notebooks/
│   ├── 01_initial_data_exploration.ipynb      # demographics, participants.tsv (Python)
│   ├── 01_initial_data_exploration_R.ipynb    # demographics, participants.tsv (R)
│   ├── 02_data_analysis.ipynb                 # early MNE exploration (Python)
│   ├── 02_data_analysis_R.ipynb               # early MNE exploration (R)
│   ├── 03_data_preparation.ipynb              # Step 2 (Python): raw EEG → per-stimulus CSVs
│   ├── 03_data_preparation_R.ipynb            # Step 2 (R): raw EEG → per-stimulus CSVs
│   ├── 04_data_preparation.ipynb              # F7 smoothing → .pkl files (Python)
│   ├── 04_data_preparation_R.ipynb            # F7 smoothing → .rds files (R)
│   ├── 05_read_data.ipynb                     # .pkl structure inspection (Python)
│   ├── 05_read_data_R.ipynb                   # .rds structure inspection (R)
│   ├── 06_plot_data.ipynb                     # visual checks of .pkl files (Python)
│   ├── 06_plot_data_R.ipynb                   # visual checks of .rds files (R)
│   ├── 07_participants_eda.ipynb              # participants EDA (Python)
│   └── 07_participants_eda_R.ipynb            # participants EDA (R)
├── Presentations/                             # Beamer slides (1st, 2nd, 3rd; compiled PDFs)
├── Practice/                                  # University lab materials and experiments
├── Makefile                                   # run `make help` for targets
├── Functional-Data-Analysis.Rproj             # RStudio project config
├── renv.lock                                  # R dependency versions
├── requirements.txt                           # Python dependencies
├── .gitignore
├── .gitmodules                                # ds006018 submodule reference
├── .renvignore                                # files/folders renv should ignore
└── README.md
```

---

## R / Python Pipeline

- R code lives in:
  - `.R` files
  - `.ipynb` notebooks that necessarily end in `_R.ipynb` (e.g.
    `03_data_preparation_R.ipynb`)
- Python code lives in:
  - `.py` files
  - `.ipynb` notebooks with regular naming (e.g. `03_data_preparation.ipynb`)

### Compatibility Notes

This project maintains two parallel pipelines — R and Python — as required. For most
analysis steps, both pipelines produce equivalent results with minor floating point
differences.

**Exception: EEG epoching step**

The per-stimuli CSVs in `ds006018_per_stimuli/` are generated by
`Notebooks/03_data_preparation.ipynb` (Python/MNE). An R equivalent
(`Notebooks/03_data_preparation_R.ipynb`, using `eegUtils`) has been implemented,
but produces amplitude values in a different scale due to a known difference between
MNE and eegUtils in how they handle BrainVision `.vhdr` calibration:

- **MNE** applies the channel resolution (0.0488281 µV/ADC unit) automatically on import
- **eegUtils** does not apply this calibration automatically

Additionally, the two libraries differ in their internal filtering and baseline
correction implementations, resulting in numerically non-identical outputs even after
manual scaling.

Both pipelines are internally consistent. All downstream R analysis (`EDA/`, `HT/`,
`REG/`) operates on the Python-generated CSVs, which are in standard µV units and
represent the ground truth for this project.

### `ds006018_functional/` and `ds006018_per_stimuli/` Coexistence

Both pipelines write their smoothed functional data objects to the same
`ds006018_functional/` directory structure, using different file extensions:

- R writes `.rds` files (`saveRDS`)
- Python writes `.pkl` files (`pickle.dump`)

The two sets of files coexist without conflict. Running either pipeline when the
other's files are already present will not overwrite or delete them.

Both pipelines also write per-stimulus CSV files to the same `ds006018_per_stimuli/`
directory structure with identical filenames and column structure
(`time, epoch, <channels>`). The intended design is that both pipelines produce
equivalent CSV outputs so that either can overwrite the other's files safely.

In practice, a known difference remains: MNE (Python) and eegUtils (R) apply
BrainVision `.vhdr` calibration differently, resulting in amplitude values that
differ in scale. All downstream analysis in both pipelines currently operates on
the Python/MNE-generated CSVs, which are in standard µV units. Resolving the R
pipeline calibration is a known open task.

### `Practice/` Folder

The `Practice/` folder contains lecture lab exercises and practice scripts used
during the course. It is not part of the main analysis pipeline and is not subject
to the dual R / Python pipeline requirement. All files there happen to be in R
(`.R` scripts and R-kernel `.ipynb` notebooks), reflecting the lecture materials
as provided.

### Known Issues

**`BandDepth` in `EDA/Smoothing_and_EDA.ipynb`**

The Band Depth (BD) plot (`23_band_depths.pdf`) shows all subjects with identical
depth values (~0.032), resulting in a flat line. Modified Band Depth (MBD) is
unaffected and matches the R output correctly. The root cause is under
investigation — likely a subtle incompatibility between scikit-fda's `BandDepth`
implementation and the specific data structure used. The MBD-based analyses
(functional boxplots, outlier detection) are not affected.

---

## Make targets

```bash
$ make help                   # Show all available targets
$ make all                    # Run the full pipeline from scratch
$ make deps                   # Install Python and R dependencies
$ make export-requirements    # Export uv dependencies to requirements.txt
$ make import-requirements    # Import requirements.txt into uv
$ make data                   # Acquire raw EEG data via datalad
$ make stimuli                # Raw EEG → per-stimulus CSVs (~40 min)
$ make functional             # Generate F7 .rds files (optional, ~40 min)
$ make assemble               # CSVs → subject matrix CSV
$ make eda                    # Smoothing + full EDA
$ make hypothesis_testing     # Run all hypothesis testing scripts
$ make regression             # Run regression analysis
$ make presentation_1         # Compile LaTeX slides for 1st presentation
$ make presentation_2         # Compile LaTeX slides for 2nd presentation
$ make presentation_3         # Compile LaTeX slides for 3rd presentation
$ make report                 # Compile final LaTeX report PDF
$ make clean                  # Remove generated outputs
$ make distclean              # Clean + remove all generated data folders
$ make clean-env              # Remove Python venv and R library
```

Full pipeline from scratch:

```bash
$ make all
```
