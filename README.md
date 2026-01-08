# 🧬✨ DNA Methylation: t-SNE & CNV Pipeline (R)

Welcome, Fahmi is here! 👋  

  
This repository is a **reproducible, helper-driven R pipeline** for analyzing **Illumina DNA methylation array (IDAT)** data, with two main superpowers:  
  
🌀 **t-SNE embedding** for contextualizing samples  
🧬 **CNV inference** directly from methylation intensities  

  
> Think of this repo as a **methylation microscope** 🔬  
> It helps you *look*, *compare*, and *contextualize* — not magically diagnose.  

---

## 🧠 What This Repo Is (and Isn’t)

### ✅ What it DOES

* Read Illumina IDAT files (450k / EPIC / EPICv2)
* Build clean beta matrices (samples × probes)
* Embed samples using **PCA → t-SNE**
* Visualize:
  * 🌍 global methylation landscapes
  * 🔍 zoomed neighborhoods
* Infer **copy number variation (CNV)** from methylation intensities
* Export **publication-ready figures & tables**
* Lock all package versions with **`renv`**
* Separate **biology decisions** from **code logic** (this is important!)

### ❌ What it very intentionally does NOT do

* RNA-seq analysis ❌
* Fusion detection ❌
* FASTQ / BAM handling ❌
* Expression inference ❌

Why?
Because **IDAT files cannot do those things** — and pretending otherwise is bad science 😌

---

## ✨ Repo Philosophy (aka “Why this feels nice to work with”)

> **The workflow is reproducible.**  
> **The helpers are reusable.**  
> **The biology is explicit.**  
> **The judgment stays with the user.**  

No hidden magic.  
No silent assumptions.  
No “trust me bro” CNV results.  

---

## 📁 Repository Structure

```text
dna-methylation-tsne-cnv/
├─ README.md
├─ renv.lock
├─ .gitignore
├─ .renvignore
├─ .Rprofile
│
├─ filter/                  # probe exclusion lists
│  ├─ amb_3965probes.vh20151030.txt
│  ├─ epicv1B2_32260probes.vh2016.txt
│  ├─ snp_7998probes.vh20151030.txt
│  └─ xy_11551probes.vh20151030.txt
│
├─ helpers/                 # reusable brain cells 🧠
│  └─ MNPTraining/          # source from another Github, usable
│  │  ├─ MNPprocessIDAT_functions.R
│  │  └─ RSpectra_pca.R
│  ├─ batch_tools.R
│  ├─ cnv_tools.R
│  ├─ config_tools.R
│  ├─ plot_tsne_tools.R
│  └─ probe_id_tools.R
│
├─ ori_script/               # original exploratory scripts (archive)
│  ├─ ori_preprocessing.R
│  ├─ ori_tsne.R
│  └─ ori_cnv.R
│
├─ src/                      # the actual pipeline 🚀
│  ├─ 01_prework.R
│  ├─ 02_embedding_tsne.R
│  └─ 03_cnv_conumee2.R
│
└─ results/                  # outputs (gitignored, always)
```

---

## 🧠 How to Think About This Repo

### 🟦 `src/` — *the workflow*

Short, readable, boring (boring = good).

* **`01_prework.R`**
  🧱 Builds the foundation

  * reads IDATs
  * preprocesses methylation
  * harmonizes probes
  * saves clean betas

* **`02_embedding_tsne.R`**
  🌌 Makes things pretty (and interpretable)

  * PCA → t-SNE
  * global & zoomed plots
  * visualization only (no biological claims)

* **`03_cnv_conumee2.R`**
  🧬 Where CNV happens

  * choose controls
  * match GEO → IDATs
  * run CNV per case
  * export everything neatly

### 🟩 `helpers/` — *the brains*

Reusable logic lives here:

* probe ID sanity
* plotting helpers
* CNV mechanics
* early “contract checks” that scream before mistakes happen 🔥

### 🟨 `ori_script/` — *the fossil record*

Old exploratory scripts, kept for transparency and provenance.

---

## 🔬 Dataset-Aware by Design

Some things **cannot** be universal in CNV analysis — and that’s okay.

| Item               | Why you must define it     |
| ------------------ | -------------------------- |
| Control classes    | Biology depends on dataset |
| GEO labels         | Vary across studies        |
| Zoom regions       | Pure visualization choice  |
| CNV interpretation | Always biological judgment |

This repo **forces those choices to be explicit**, instead of hiding them.

> Reproducible ≠ pretending biology is generic.

---

## 🔒 Data Policy & Privacy

* 🔐 **Raw IDAT files are NEVER committed**
* 📁 `results/` is always gitignored
* 🧪 Only code + dummy metadata live here
* 🚫 No patient identifiers, ever
  
This repo is safe to share **publicly as code only**.

---

## 🔁 Reproducibility with `renv`

```r
install.packages("renv")
renv::restore()
```

That’s it.  
Everyone now runs the **same R universe** 🌍  

---

## ▶️ How to Run

```r
source("src/01_prework.R")
source("src/02_embedding_tsne.R")
source("src/03_cnv_conumee2.R")
```

Grab coffee ☕  
Check `results/` 📂  
Smile 😄  

---

## 📤 What You Get Out

### 🌀 t-SNE

* Global methylation landscape
* Zoomed neighborhood plots
* Sample-level coordinates

### 🧬 CNV

* Genome-wide CNV plots
* Segment tables
* Focal CNV exploration (EGFR, CDKN2A/B, etc.)
* Clean per-case output bundles

---

## 🧪 Intended Use

This repository is for:

* Exploratory methylation-based tumor analysis
* Case contextualization vs public references
* CNV screening from methylation arrays
* Research & education

It is **not** a diagnostic tool.

---

## 🙏 Acknowledgements & Upstream References

This repository builds upon and adapts ideas, helper functions, and filter resources from the following upstream project:

**MNP Training Repository**  
🔗 [https://github.com/mwsill/mnp_training](https://github.com/mwsill/mnp_training)

Specifically, inspiration and adapted logic were taken from:

- `MNPprocessIDAT_functions.R`  
  → preprocessing helpers for Illumina methylation IDAT files ==> In our repo, we add some modification for EPICv2 
- `RSpectra_pca.R`  
  → fast PCA routines for large methylation matrices  
- `filter/` resources  
  → probe exclusion lists (SNP, sex chromosomes, ambiguous probes)

These components have been **integrated, refactored, and extended** to support:
- EPICv2 compatibility
- modern `renv`-based reproducibility
- modular helper-based workflows
- explicit dataset-aware configuration

All mistakes or modifications in this repository are the responsibility of the author.

---

## ✍️ Author

**Dr. Muhammad Nur Fahmi**
