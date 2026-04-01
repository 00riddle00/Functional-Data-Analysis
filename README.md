<!-- vim: set ft=markdown fenc=utf-8 tw=88 nu ai si et ts=2 sw=2: -->

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
- [Notes](#notes)
- [Make targets](#make-targets)

## Quick start

```bash
git clone --recurse-submodules https://github.com/OpenNeuroDatasets/ds006018
make deps
make all
```

For individual steps, run `make help`.

---

## Prerequisites

### Git

- **Linux:** `sudo apt install git` (Ubuntu/Debian) or use your distro's package manager
- **macOS:** `brew install git`
- **Windows:** https://git-scm.com/download/win — after install, run: `git config
  --global core.longpaths true`

### Python (3.14+)

- **Linux:** `sudo apt install python3 python3-venv python3-pip` (or use your distro's
  package manager)
- **macOS:** `brew install python`
- **Windows:** https://www.python.org — check "Add to PATH" during install

### R (4.5+)

- **Linux:** https://cran.r-project.org/bin/linux/ (follow distro-specific instructions)
- **macOS:** https://cran.r-project.org/bin/macosx/
- **Windows:** https://cran.r-project.org/bin/windows/base/

### DataLad and git-annex

Needed to acquire the raw 10 GB EEG data.

- **Linux:** `sudo apt install datalad git-annex`
- **macOS:** `brew install git-annex && pip install datalad`
- **Windows:** Install git-annex from https://git-annex.branchable.com/install/Windows/
  then `pip install datalad`.

### LaTeX

Needed to compile the presentation locally.

- **Linux:** `sudo apt install texlive-full` (or minimal:
  `texlive-base texlive-latex-recommended texlive-latex-extra
  texlive-fonts-recommended texlive-xetex latexmk`)
- **macOS:** https://www.tug.org/mactex/ or `brew install --cask mactex`
- **Windows:** https://miktex.org/download (select "Install missing packages
  on the fly"). Alternatively, TeX Live: https://www.tug.org/texlive/windows.html

### System libraries (Linux only)

Some R packages require system libraries. Install before running `renv::restore()`:

```bash
sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
  libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```

### Jupyter R kernel (optional)

Only needed to run R-based notebooks in Jupyter.

```r
install.packages("IRkernel")
IRkernel::installspec(
  name = paste0("r-", basename(getwd())),
  displayname = paste0("R (", basename(getwd()), ")")
)
```

This creates a project-specific kernel name. Renv activation is handled
by the first cell in each R notebook.

### Make (optional)

Needed if you want to use the Makefile for automation.

- **Linux:** Usually pre-installed. If not: `sudo apt install make`
- **macOS:** Included with Xcode command line tools (`xcode-select --install`) or just
  `brew install make`
- **Windows:** Comes with Rtools4 (which you likely have for R package compilation).
  Otherwise: `choco install make`. Run from Git Bash or Rtools terminal, not PowerShell.

---

## Setup

### Clone

```bash
git clone --recurse-submodules https://github.com/OpenNeuroDatasets/ds006018
```

If already cloned without `--recurse-submodules`:

```bash
git submodule update --init
```

### Install dependencies

**Python:**

```bash
python3 -m venv .venv
source .venv/bin/activate        # Linux/macOS
# .venv\Scripts\activate         # Windows
pip install -r requirements.txt
```

**R:**

```r
install.packages("renv")
renv::restore()
```

---

## Pipeline

The pipeline has 4 steps. Each depends on the output of the previous one.
Run `make all` to execute steps 3–4 and compile the presentation, or run
each step individually.

### Step 1: Acquire raw EEG data (10+ GB, takes a while)

**DataLad:**

```bash
cd ds006018 && datalad get . && cd ..
```

### Step 2: Generate per-stimulus CSVs (~40 min)

```bash
cd notebooks
jupyter nbconvert --to notebook --execute 03_data_preparation.ipynb --output /dev/null
cd ..
```

Reads raw BrainVision EEG files, filters, epochs, and exports one CSV per subject per
task per stimulus.

- **Input:** `ds006018/sub-XXX/eeg/*.vhdr`
- **Output:** `ds006018_per_stimuli/sub-XXX/task-XXX_Stimulus_SXX.csv` (50+ GB)

### Step 2b (optional): Generate smoothed .rds files (~40 min)

Only needed if you want F7-channel functional data objects. The main analysis (Steps
3–4) does not use these — it extracts FC1 directly from the CSVs.

```bash
cd Notebooks
jupyter nbconvert --to notebook --execute 04_data_preparation_R.ipynb --output /dev/null
cd ..
```

- **Input:** `ds006018_per_stimuli/sub-XXX/*.csv`
- **Output:** `ds006018_functional/sub-XXX/*.rds` (45+ MB, F7 channel only)

**Note:** Some output filenames contain `:` characters (e.g., `LostSamples:264.rds`)
which are invalid on Windows. The last cell of the notebook renames these automatically.

### Step 3: Assemble the across-subjects data matrix

```bash
cd EDA && Rscript assemble_subject_flanker_S2_FC1.R && cd ..
```

Extracts FC1 channel from flanker S2 stimulus, averages epochs if needed,
combines all 62 subjects into one CSV.

- **Input:** `ds006018_per_stimuli/sub-XXX/task-flanker_Stimulus_S2.csv`
- **Output:** `EDA/Flanker_stimulus_FC1_channel.csv` + `EDA/subject_metadata.csv`

### Step 4: Smoothing + Exploratory Data Analysis (EDA)

```bash
cd EDA && Rscript Smoothing_and_EDA.R && cd ..
```

Applies B-spline smoothing with GCV-selected λ, then runs the full functional
EDA: mean/SD, covariance, FPCA, depth, outlier detection, boxplots, rainbow
plots.

- **Input:** `EDA/Flanker_stimulus_FC1_channel.csv`
- **Output:** 30+ PDF plots + 1 text file + `fd_smooth.rds` in `EDA/outputs/`

### Step 5 (optional): Compile presentation

```bash
cp EDA/outputs/*.pdf presentation/
cd presentation && latexmk -xelatex main.tex && cd ..
```

The presentation uses Beamer with XeLaTeX. If `latexmk` is unavailable, run
`xelatex main.tex` twice manually.

### LaTeX troubleshooting

- **xeCJK error:** If you don't need CJK (Chinese/Japanese/Korean) support,
  remove `\usepackage{xeCJK}` from `main.tex`.
- **Font error on Linux:** Install Utopia font: `sudo apt install
  texlive-fonts-extra`, or remove `\usepackage{utopia}` from `main.tex`.
- **Missing packages on Windows (MiKTeX):** MiKTeX installs missing packages
  automatically on first compile. If prompted, click "Install".
- **Missing packages on Linux (TeX Live):** `tlmgr install <package-name>`

---

## Folder structure

```
.
├── ds006018/                      # Git submodule — raw EEG data (10+ GB via datalad)
├── ds006018_per_stimuli/          # .gitignore — intermediate CSVs (50+ GB)
├── ds006018_functional/           # .gitignore — F7 .rds files (unused)
├── EDA/
│   ├── assemble_subject_flanker_S2_FC1.R  # Step 3: CSVs → subject matrix
│   ├── Smoothing_and_EDA.R        # Step 4: smoothing + full EDA
│   ├── Flanker_stimulus_FC1_channel.csv  # 62 subjects × 501 time points
│   ├── subject_metadata.csv       # epoch counts per subject
│   └── outputs/                   # All PDFs, fd_smooth.rds, text results
├── Notebooks/
│   ├── 01_initial_data_exploration.ipynb  # demographics, participants.tsv
│   ├── 02_data_analysis.ipynb             # early MNE exploration
│   ├── 03_data_preparation.ipynb          # Step 2: raw EEG → per-stimulus CSVs
│   ├── 04_data_preparation_R.ipynb        # F7 smoothing → .rds files
│   ├── 05_read_data_R.ipynb               # .rds structure inspection
│   └── 06_plot_data_R.ipynb               # visual checks of .rds files
├── Presentations/                 # Beamer slides (main.tex, compiled PDFs)
├── Slides/                        # University lecture slides (reference only)
├── Practice/                      # University lab materials and our experiments
├── Makefile                       # run `make help` for targets
├── Functional-Data-Analysis.Rproj # RStudio project config
├── renv.lock                      # R dependency versions
├── requirements.txt               # Python dependencies
├── .gitignore
├── .gitmodules                    # ds006018 submodule reference
├── .renvignore                    # files/folders renv should ignore
└── README.md
```

---

## Notes

- `ds006018/` is a git submodule pointing to
  https://github.com/OpenNeuroDatasets/ds006018.git. Never commit large files into it.
- `ds006018_per_stimuli/` is in `.gitignore` (54 GB). Regenerate from step 2.
- `ds006018_functional/` is in `.gitignore`. Contains F7-channel .rds files. Currently
  unused — the main analysis uses FC1 extracted directly from the CSVs.
- The analysis uses **channel FC1** (frontal-central, relevant for attentional
  conflict processing in the Flanker task) and **stimulus S2**.
- 62 subjects had flanker S2 data: 55 with 1 trial, 7 with 2 trials (averaged).

---

## Make targets

```bash
make help          # show all available targets
make deps          # install Python and R dependencies
make data          # acquire raw EEG data via datalad
make stimuli       # raw EEG → per-stimulus CSVs (~40 min)
make assemble      # CSVs → subject matrix CSV (~30 sec)
make eda           # smoothing + full EDA (~2 min)
make presentation  # compile LaTeX slides
make all           # run eda + presentation
make clean         # remove generated outputs
make clean-all     # also remove per-stimulus CSVs
```

Full pipeline from scratch:

```bash
make deps data stimuli assemble eda presentation
```

