# Protein Structure Analysis

`protein_structure_analysis.py` analyzes a PDB structure and reports:

- **Primary structure**: chain sequences, amino acid composition, estimated molecular weight.
- **Secondary structure**: helix/sheet/coil fractions (from phi/psi heuristic).
- **Tertiary structure**: radius of gyration, contact count, residue counts.
- **Quaternary structure**: inferred oligomeric state and inter-chain contacts.
- **Binding pocket analysis**: cavity candidates from a grid-based heuristic.
- **Tunnel prediction**: approximate channel directions and bottlenecks from pocket geometry.

## Install

```bash
pip install biopython numpy scipy
```

## Usage

```bash
python protein_structure_analysis.py --pdb input.pdb --out analysis.json
```

## Notes

- Binding pocket and tunnel methods are lightweight geometric heuristics.
- For production-grade cavity/tunnel characterization, use dedicated tools (e.g., fpocket, CAVER) and compare against this script.


## Online-style UI (shareable link)

You can run a browser UI and get a temporary public link:

```bash
pip install gradio biopython numpy scipy
python app.py
```

When the app starts, Gradio prints two URLs:
- a local URL (for your machine)
- a temporary public `https://*.gradio.live` URL that you can share.


## Deploy as a permanent online tool (Render)

This repo now includes `requirements.txt`, `Procfile`, and `render.yaml` for deployment.

### Option A: Render Blueprint (recommended)
1. Push this repo to GitHub.
2. In Render, choose **New +** → **Blueprint**.
3. Select your repo; Render detects `render.yaml`.
4. Click **Apply** and wait for build.
5. Render gives you a permanent URL like `https://protein-structure-analyzer.onrender.com`.

### Option B: Generic Python host
Use start command:

```bash
python app.py --server-name 0.0.0.0 --server-port $PORT --no-share
```

`--no-share` disables temporary tunnels and is better for hosted deployments.
