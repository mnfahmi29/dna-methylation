**⚠️⚠️Desclaimer: due to my new step in Github, I will post my code as soon as possible⚠️⚠️**
# DNA Methylation t-SNE & CNV Analysis (R)

A reproducible R pipeline for **DNA methylation array (IDAT)** analysis focused on:

- **t-SNE embedding** against a public reference cohort for tumor classification / contextualization  
- **Copy Number Variation (CNV)** analysis derived from methylation intensities

This repository is intentionally scoped to **methylation-based analyses only**.  
It does **not** perform RNA-seq, fusion detection, or transcriptomic analyses.

---

## ✨ Features

- Read and preprocess Illumina methylation IDAT files
- Harmonize case samples with a public reference cohort
- Perform dimensionality reduction (PCA → t-SNE)
- Generate:
  - Global t-SNE embedding
  - Zoomed-in neighborhood plots for specific tumor classes
- Perform CNV analysis using methylation intensities
- Export publication-ready figures and summary tables
- Fully reproducible via `renv`

---

## 📁 Repository Structure

```text
dna-methylation-tsne-cnv-r/
├─ README.md
├─ renv.lock
├─ .gitignore
├─ R/
│  ├─ 01_setup.R
│  ├─ 02_idat_import_controls.R
│  ├─ 03_tsne_classification.R
│  ├─ 04_cnv_conumee2.R
│  ├─ 05_reports_export.R
│  └─ utils/
│     ├─ io_helpers.R
│     ├─ plotting_helpers.R
│     └─ safe_filename.R
├─ data/
│  ├─ raw/            # local only (IDATs not committed)
│  └─ example/        # dummy metadata only
├─ results/           # generated outputs (ignored by git)
└─ docs/
   ├─ figures/
   └─ notes.md
````

---

## 🔒 Data Policy & Privacy

* **Raw IDAT files are never committed**
* `data/raw/` is intentionally gitignored
* Only **dummy / synthetic metadata** is stored in `data/example/`
* No patient identifiers (name, MRN, DOB, hospital number) should appear anywhere in the repo

This repository is safe to share publicly **as code only**.

---

## 🔁 Reproducibility (`renv`)

This project uses **`renv`** to lock R package versions.

### First-time setup

```r
install.packages("renv")
renv::restore()
```

This installs the exact package versions used to develop the pipeline.

---

## ▶️ How to Run

After placing your IDAT files locally and preparing a sample sheet:

```r
source("R/01_setup.R")
source("R/02_idat_import_controls.R")
source("R/03_tsne_classification.R")
source("R/04_cnv_conumee2.R")
source("R/05_reports_export.R")
```

All outputs will be written to the `results/` directory.

---

## 📤 Outputs

### t-SNE

* `results/tsne_global.pdf`
* `results/tsne_zoom_<CLASS>.pdf`
* Sample-level embedding coordinates

### CNV

* Genome-wide CNV segment plots
* Gene-level focal CNV plots (e.g. EGFR, CDKN2A/B, MDM2)
* Tabular summaries of focal events

---

## 🚫 Explicit Non-Goals

This repository **does not**:

* Perform RNA-seq analysis
* Detect gene fusions
* Analyze FASTQ / BAM files
* Infer transcriptomic expression

IDAT files **cannot** support these analyses, and they are intentionally excluded.

---

## 📚 Methods (high-level)

* Methylation preprocessing: Illumina array preprocessing via `minfi`
* Dimensionality reduction: PCA followed by t-SNE
* CNV inference: methylation intensity-based CNV calling
* Visualization: `ggplot2`, base plotting, and PDF export

Detailed method notes are available in `docs/`.

---

## 🧪 Intended Use

This repository is designed for:

* Exploratory methylation-based tumor classification
* Case contextualization against known reference cohorts
* CNV screening from methylation data
* Research and educational use

It is **not** a diagnostic tool.

---

## 📜 License

Specify your license here (e.g. MIT, BSD-3, or institutional license).

---

## ✍️ Author

Dr. Muhammad Nur Fahmi
