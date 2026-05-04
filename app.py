#!/usr/bin/env python3
"""Gradio app for protein structure analysis.

Run:
  python app.py

This starts a local UI and (with share=True) also prints a temporary public link.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import gradio as gr

from protein_structure_analysis import run_analysis


def analyze_uploaded_pdb(file_obj):
    if file_obj is None:
        return "Please upload a PDB file.", None

    src = Path(file_obj.name)
    if src.suffix.lower() != ".pdb":
        return "Please upload a .pdb file.", None

    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_copy = Path(tmpdir) / src.name
        pdb_copy.write_bytes(src.read_bytes())

        result = run_analysis(str(pdb_copy))
        json_text = json.dumps(result, indent=2)

        out_path = Path(tmpdir) / "analysis.json"
        out_path.write_text(json_text)
        return json_text, str(out_path)


with gr.Blocks(title="Protein Structure Analyzer") as demo:
    gr.Markdown("""
    # Protein Structure Analyzer
    Upload a **PDB** file to run:
    - Primary / Secondary / Tertiary / Quaternary analysis
    - Binding pocket analysis
    - Tunnel prediction
    """)

    pdb_file = gr.File(label="Upload PDB file", file_types=[".pdb"])
    run_btn = gr.Button("Analyze")

    json_output = gr.Code(label="Analysis JSON", language="json")
    json_download = gr.File(label="Download analysis.json")

    run_btn.click(fn=analyze_uploaded_pdb, inputs=[pdb_file], outputs=[json_output, json_download])


if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("--server-name", default="127.0.0.1")
    parser.add_argument("--server-port", type=int, default=int(os.getenv("PORT", "7860")))
    parser.add_argument("--no-share", action="store_true", help="Disable temporary public Gradio link")
    args = parser.parse_args()

    demo.launch(
        server_name=args.server_name,
        server_port=args.server_port,
        share=not args.no_share,
    )
