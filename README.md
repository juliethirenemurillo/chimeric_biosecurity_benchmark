# Chimaera Watch
### AlphaFold2-based chimeric protein viability screening

Part of **BioBouncer** — a biosecurity screening stack built on top of SecureDNA and the IBBIS Common Mechanism.

---

## What this module does

Current DNA screening catches dangerous sequences by similarity to known threats. A motivated bad actor can bypass this by ordering two individually harmless sequences and fusing them in the lab.

**Chimaera Watch** closes this gap:

```
DNA1 → protein A → sequence screen PASS → structure screen PASS ✓
DNA2 → protein B → sequence screen PASS → structure screen PASS ✓

                    ↓  (same order / same customer)

              Predict A+B fusion with AlphaFold2
                    ↓
         Is the chimera structurally viable?
                    ↓
              YES → feed into structure
                    similarity search
                    against harmful protein DB
                    ↓
              FLAG or PASS
```

This module implements the viability prediction step. It uses AlphaFold2 pLDDT scores to detect whether a chimeric fusion is structurally coherent — and if so, hands the predicted structure to the structural similarity layer for harmful-protein comparison.

---

## Where it sits in BioBouncer

```
┌─────────────────────────────────────────────────────┐
│                  BioBouncer pipeline                 │
├─────────────────────────────────────────────────────┤
│  Layer 1  │  DNA sequence screening                 │
│           │  SecureDNA + IBBIS Common Mechanism     │
├─────────────────────────────────────────────────────┤
│  Layer 2  │  Structural similarity screening        │
│           │  Foldseek + US-align vs harmful DB      │
│           │  (BioBouncer structure-screen module)   │
├─────────────────────────────────────────────────────┤
│  Layer 3  │  Chimaera Watch  ← THIS REPO            │
│           │  AlphaFold2 viability prediction        │
│           │  for multi-sequence orders              │
│           │  → feeds viable chimeras back to Layer 2│
└─────────────────────────────────────────────────────┘
```

---

## Scientific approach

For each chimera we run three AlphaFold2 predictions:

1. Domain A alone → baseline pLDDT
2. Domain B alone → baseline pLDDT
3. Full chimera A+B → pLDDT per domain in fusion context

**Key metric — delta pLDDT:**
```
Δ pLDDT = domain_pLDDT_in_chimera − domain_pLDDT_alone

Large negative Δ → fusion disrupts folding → chimera NOT viable → no threat
Near-zero Δ      → domains fold independently → chimera viable → escalate to Layer 2
```

pLDDT is AlphaFold2's per-residue confidence score (0–100). It is read directly from the B-factor column of the predicted PDB structure — no additional tools required.

---

## Test cases

| ID | Chimera | Expected | Basis |
|----|---------|----------|-------|
| E001 | EGFP + mCherry | WORKS | Campbell et al. 2002 — known tandem fusion |
| F001 | Hsp70 (Ssa1) + GFP | FAILS | Published experimental failure |
| F002 | IHFbeta + GFP (343 aa) | FAILS | Smallest case; IHFbeta is a dimer in nature |
| F003 | p53(Y220C) + GFP | FAILS | Published experimental failure |

**Benchmark logic:**
- E001 should show small Δ pLDDT (domains preserved in fusion)
- F001–F003 should show large negative Δ pLDDT (fusion disrupts folding)

If AlphaFold2 correctly distinguishes these, the benchmark validates that structural viability prediction can detect chimeric evasion attempts.

---

## Output

The module produces a JSON verdict per chimera order:

```json
{
  "order_ids": ["order_001", "order_002"],
  "sequences": ["IHFbeta", "GFP"],
  "length_aa": 343,
  "individual_pass": true,
  "domain_scores": [
    {
      "domain": "IHFbeta",
      "residues": "1-99",
      "pLDDT_alone": 83.0,
      "pLDDT_chimera": 55.0,
      "delta": -28.0
    },
    {
      "domain": "GFP",
      "residues": "106-343",
      "pLDDT_alone": 96.0,
      "pLDDT_chimera": 58.0,
      "delta": -38.0
    }
  ],
  "chimera_verdict": "FLAG",
  "threshold_used": -10,
  "action": "send_to_structural_similarity_search"
}
```

When `chimera_verdict` is `FLAG`, the predicted chimera PDB is passed to the BioBouncer structural similarity layer (Layer 2) for comparison against the harmful protein database.

---

## Results so far

| Sequence | Mean pLDDT | Source |
|----------|-----------|--------|
| GFP alone | **96.04** | AF2 measured ✅ |
| EGFP alone | ~95 | Running |
| mCherry alone | ~88 | Running |
| IHFbeta alone | ~83 | Running |
| E001 chimera | pending | Running |
| F002 chimera | pending | Running |

GFP alone result confirms AF2 is predicting correctly — 96.2% of residues above pLDDT 70, all 5 models consistent.

---

## Quickstart

```bash
git clone https://github.com/juliethirenemurillo/chimeric_biosecurity_benchmark
cd chimeric_biosecurity_benchmark
git checkout bio/test-cases

# Install dependencies
pip install -r requirements.txt

# Run pLDDT extraction on completed AF2 results
python scripts/extract_plddt.py

# Run full benchmark scoring (uses placeholders for pending jobs)
python scripts/benchmark_scoring.py
```

### Running AlphaFold2 on CSD3 (Cambridge HPC)

```bash
# Submit baseline jobs
bash scripts/submit_af2_jobs.sh --dry-run   # preview
bash scripts/submit_af2_jobs.sh             # submit

# For sequences with known CIF issues on CSD3
sbatch --export=ALL,JOB_NAME=...,FASTA_PATH=...,OUTPUT_DIR=... \
  scripts/af2_retry.sh
```

**Requirements:** AlphaFold2 2.3.2, icelake partition, BRYANT-SL3-CPU account

---

## Repository structure

```
chimeric_biosecurity_benchmark/
├── chimeras/                  # FASTA files for all test cases
│   ├── GFP_alone.fasta
│   ├── EGFP_alone.fasta
│   ├── mCherry_alone.fasta
│   ├── IHFbeta_alone.fasta
│   ├── E001.fasta             # EGFP+mCherry (WORKS)
│   ├── F002.fasta             # IHFbeta+GFP  (FAILS)
│   └── ...
├── scripts/
│   ├── af2_single.sh          # SLURM job script (single AF2 run)
│   ├── af2_retry.sh           # SLURM job script with CIF retry loop
│   ├── submit_af2_jobs.sh     # Launcher for parallel submissions
│   ├── extract_plddt.py       # Per-residue pLDDT extraction + plots
│   └── benchmark_scoring.py   # Full benchmark table + figures
├── af2_results/               # AF2 output (gitignored, on HPC)
├── results/                   # Plots and CSVs
└── run_pipeline.py            # End-to-end pipeline runner
```

---

## Team

Part of the **BioBouncer** project — LISA Hackathon 2025.

- **Chimaera Watch (this repo):** Julieth Murillo — chimeric viability prediction via AlphaFold2
- **Structure screening:** [colleague] — Foldseek structural similarity vs harmful protein DB
- **Pipeline integration:** [developers] — IBBIS / benchtop synthesiser integration
- **Policy:** [policy colleague] — biosecurity governance and operator verification

Related repo: [BioBouncer structure-screen](https://github.com/[colleague-repo])

---

## Biosecurity relevance

If AlphaFold2 correctly predicts chimeric viability:
- A screener can flag structurally coherent multi-part orders before synthesis
- Current sequence-only screening would miss these combinations entirely
- This supports the case for **function-based screening** as a complement to sequence similarity

> Sequence similarity is not sufficient for screening in the age of AI-generated biology.
