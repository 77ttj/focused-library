from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Callable, Iterable

import gradio as gr
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

import generation_by_compared_pocket_similarity as compared_pocket_generator
import generation_by_keyaas_overlap_similarity_specify_keyaas as keyaas_generator
import generation_by_pocket_similarity as pocket_generator
import generation_by_seq_similarity_sequence_only as sequence_generator


NUM_REFERENCE_POCKETS = 1000
MAX_GALLERY_MOLECULES = 4


CSS = """
.gradio-container {
    width: min(100%, 980px) !important;
    max-width: 980px !important;
    margin: 0 auto !important;
    padding: 12px 18px 28px !important;
    background: #ffffff !important;
}
#app-title {
    text-align: center;
    margin: 0 0 14px 0;
}
#app-title h1 {
    margin: 0;
    font-size: 26px;
    line-height: 1.25;
    font-weight: 700;
    color: #111111;
}
#method-nav {
    display: flex !important;
    flex-wrap: nowrap !important;
    gap: 8px !important;
    margin: 0 0 16px 0 !important;
}
.method-tab {
    flex: 1 1 0 !important;
    min-width: 0 !important;
    min-height: 48px !important;
    padding: 9px 10px !important;
    border: 1px solid #eadbe9 !important;
    border-radius: 8px !important;
    background: #f7eef8 !important;
    color: #111111 !important;
    font-size: 15px !important;
    font-weight: 600 !important;
    box-shadow: none !important;
    white-space: nowrap !important;
}
.method-tab.selected,
.method-tab.primary {
    border-color: #e5d3e6 !important;
    background: #f0e4f1 !important;
    color: #111111 !important;
}
.method-tab:hover {
    background: #eadbea !important;
}
.method-panel {
    border: 0 !important;
    padding: 0 !important;
}
.method-panel .form {
    background: transparent !important;
}
.input-row {
    display: flex !important;
    flex-wrap: nowrap !important;
    align-items: stretch !important;
    gap: 16px !important;
}
.input-fields,
.input-action {
    min-width: 0 !important;
}
.input-action {
    display: flex !important;
    align-items: flex-start !important;
}
.action-button {
    width: 100% !important;
    min-height: 86px !important;
    border: 1px solid #eadbe9 !important;
    border-radius: 8px !important;
    background: #f0e4f1 !important;
    color: #111111 !important;
    font-size: 20px !important;
    font-weight: 700 !important;
    box-shadow: none !important;
}
.action-button:hover {
    background: #eadbea !important;
}
.save-row {
    display: flex !important;
    flex-wrap: nowrap !important;
    gap: 16px !important;
}
.save-button {
    min-height: 66px !important;
    border: 1px solid #eadbe9 !important;
    border-radius: 8px !important;
    background: #f0e4f1 !important;
    color: #111111 !important;
    font-size: 20px !important;
    font-weight: 700 !important;
    box-shadow: none !important;
}
.save-button:hover {
    background: #eadbea !important;
}
.method-panel [data-testid="block-info"] {
    padding: 0 0 7px 0 !important;
    border-radius: 0 !important;
    background: transparent !important;
    color: #6b7280 !important;
    font-size: 13px !important;
    font-weight: 400 !important;
}
.method-panel input,
.method-panel textarea,
.method-panel .input-container {
    border-color: #d1d5db !important;
    box-shadow: none !important;
}
.method-panel input:focus,
.method-panel textarea:focus {
    border-color: #c4a7c8 !important;
}
.output-textbox textarea {
    font-family: Consolas, "Courier New", monospace !important;
    font-size: 12px !important;
}
.molecule-gallery {
    border: 0 !important;
    background: transparent !important;
}
.molecule-gallery .grid-wrap {
    background: #ffffff !important;
}
footer {
    display: none !important;
}
@media (max-width: 560px) {
    #method-nav {
        flex-wrap: wrap !important;
    }
    .method-tab {
        flex: 1 1 calc(50% - 4px) !important;
    }
    .input-row,
    .save-row {
        flex-wrap: wrap !important;
    }
}
"""


def _parse_key_aas(value: str) -> list[str]:
    residues = [item.upper() for item in re.split(r"[\s,;]+", value.strip()) if item]
    if not residues:
        raise gr.Error("Please enter at least one key amino acid.")
    invalid = [item for item in residues if not re.fullmatch(r"[A-Z]{3}", item)]
    if invalid:
        raise gr.Error(f"Invalid amino acid code(s): {', '.join(invalid)}")
    return residues


def _filter_by_minimum_size(smiles_list: Iterable[str], minimum_size: float) -> list[str]:
    filtered = []
    seen = set()
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        if Descriptors.MolWt(mol) < float(minimum_size):
            continue
        canonical = Chem.MolToSmiles(mol, canonical=True)
        if canonical not in seen:
            seen.add(canonical)
            filtered.append(canonical)
    return filtered


def _draw_molecules(smiles_list: list[str]) -> list:
    images = []
    for smiles in smiles_list[:MAX_GALLERY_MOLECULES]:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            images.append(Draw.MolToImage(mol, size=(320, 220)))
    return images


def _run_generator(
    mode: str,
    primary_input,
    minimum_size: float,
):
    if minimum_size is None:
        minimum_size = 0

    if mode == "sequence":
        sequence = str(primary_input or "").strip().upper()
        if not sequence:
            raise gr.Error("Please enter a pocket sequence.")
        runner: Callable = lambda: sequence_generator.process(sequence, NUM_REFERENCE_POCKETS)
    elif mode == "keyaas":
        key_aas = _parse_key_aas(str(primary_input or ""))
        runner = lambda: keyaas_generator.process(key_aas, NUM_REFERENCE_POCKETS)
    elif mode == "pocket":
        source_pdb = str(primary_input or "").strip()
        if not source_pdb:
            raise gr.Error("Please upload a target PDB file.")
        if not Path(source_pdb).is_file():
            raise gr.Error("The selected PDB file does not exist.")
        runner = lambda: pocket_generator.process(source_pdb, NUM_REFERENCE_POCKETS)
    elif mode == "paired_pocket":
        source_pdb = str(primary_input or "").strip()
        if not source_pdb:
            raise gr.Error("Please upload a target PDB file.")
        if not Path(source_pdb).is_file():
            raise gr.Error("The selected PDB file does not exist.")
        runner = lambda: compared_pocket_generator.process(source_pdb, NUM_REFERENCE_POCKETS)
    else:
        raise gr.Error(f"Unknown generator: {mode}")

    try:
        generated = runner()
    except gr.Error:
        raise
    except Exception as exc:
        raise gr.Error(f"Generation failed: {exc}") from exc

    filtered = _filter_by_minimum_size(generated or [], minimum_size)
    if not filtered:
        raise gr.Error("No molecules passed the current minimum molecule size.")

    gr.Info(f"Generated {len(filtered)} molecules.")
    return filtered, repr(filtered[:MAX_GALLERY_MOLECULES]), _draw_molecules(filtered)


def _save_smiles(smiles_list: list[str], save_path: str):
    if not smiles_list:
        raise gr.Error("Generate molecules before saving.")
    if not save_path or not save_path.strip():
        raise gr.Error("Please enter a save path.")

    output_path = Path(save_path.strip()).expanduser()
    if not output_path.is_absolute():
        output_path = Path(__file__).resolve().parent / output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(smiles_list) + "\n", encoding="utf-8")
    gr.Info(f"Saved {len(smiles_list)} SMILES to {output_path}")
    return gr.update(value=str(output_path))


def _build_output_panel(prefix: str):
    generated_state = gr.State([])
    generated_smiles = gr.Textbox(
        label="Generated SMILES",
        lines=4,
        interactive=False,
        elem_classes=["output-textbox"],
    )
    molecule_gallery = gr.Gallery(
        show_label=False,
        columns=MAX_GALLERY_MOLECULES,
        rows=1,
        height=225,
        object_fit="contain",
        preview=False,
        elem_classes=["molecule-gallery"],
    )
    with gr.Row(elem_classes=["save-row"]):
        save_path = gr.Textbox(
            label="Save Path",
            value="output.smi",
            lines=1,
            max_lines=1,
            scale=2,
        )
        save_button = gr.Button(
            "Save SMILES",
            elem_classes=["save-button"],
            scale=1,
        )

    save_button.click(
        fn=_save_smiles,
        inputs=[generated_state, save_path],
        outputs=[save_path],
        api_name=f"save_smiles_{prefix}",
    )
    return generated_state, generated_smiles, molecule_gallery


def _build_generation_tab(
    mode: str,
    prefix: str,
    primary_input_factory: Callable,
):
    with gr.Row(equal_height=True, elem_classes=["input-row"]):
        with gr.Column(scale=2, elem_classes=["input-fields"]):
            primary_input = primary_input_factory()
            minimum_size = gr.Slider(
                label="Minimum Molecule Size",
                minimum=0,
                maximum=500,
                value=200,
                step=1,
            )
        with gr.Column(scale=1, elem_classes=["input-action"]):
            run_button = gr.Button(
                "Recap and Combination",
                elem_classes=["action-button"],
            )

    generated_state, generated_smiles, molecule_gallery = _build_output_panel(prefix)
    run_button.click(
        fn=lambda value, minimum: _run_generator(mode, value, minimum),
        inputs=[primary_input, minimum_size],
        outputs=[generated_state, generated_smiles, molecule_gallery],
        api_name=f"generate_{prefix}",
        concurrency_limit=1,
        concurrency_id="focused_library_generation",
        show_progress="minimal",
    )
    return primary_input, minimum_size


def _switch_method(selected_method: str, panels, buttons):
    panel_updates = [
        gr.update(visible=name == selected_method)
        for name, _panel in panels
    ]
    button_updates = [
        gr.update(
            elem_classes=[
                "method-tab",
                *(["selected"] if name == selected_method else []),
            ],
            variant="primary" if name == selected_method else "secondary",
        )
        for name, _button in buttons
    ]
    return [*panel_updates, *button_updates]


def _build_method_panel(
    mode: str,
    prefix: str,
    primary_input_factory: Callable,
    visible: bool,
):
    with gr.Column(
        visible=visible,
        elem_id=f"panel-{prefix}",
        elem_classes=["method-panel"],
    ) as panel:
        _build_generation_tab(mode, prefix, primary_input_factory)
    return panel


def create_demo() -> gr.Blocks:
    with gr.Blocks(
        title="Focused-Library Generator",
        css=CSS,
        theme=gr.themes.Base(),
    ) as demo:
        gr.HTML(
            '<div id="app-title"><h1>🧪 Focused-Library Generator 🌿</h1></div>'
        )

        method_specs = [
            ("sequence", "Sequence Similarity"),
            ("keyaas", "KeyAA Similarity"),
            ("pocket", "Pocket Similarity"),
            ("paired_pocket", "Paired-pocket Similarity"),
        ]

        with gr.Row(elem_id="method-nav"):
            nav_buttons = []
            for index, (name, label) in enumerate(method_specs):
                nav_buttons.append(
                    (
                        name,
                        gr.Button(
                            label,
                            variant="primary" if index == 0 else "secondary",
                            elem_classes=[
                                "method-tab",
                                *(["selected"] if index == 0 else []),
                            ],
                        ),
                    )
                )

        panels = [
            (
                "sequence",
                _build_method_panel(
                    "sequence",
                    "sequence",
                    lambda: gr.Textbox(
                        label="Sequence Input",
                        value="FYSYLGAMFFLL",
                        lines=3,
                    ),
                    visible=True,
                ),
            ),
            (
                "keyaas",
                _build_method_panel(
                    "keyaas",
                    "keyaas",
                    lambda: gr.Textbox(
                        label="Key Amino Acids",
                        value="TYR, TYR, TYR",
                        placeholder="Example: TYR, TYR, TYR",
                        lines=3,
                    ),
                    visible=False,
                ),
            ),
            (
                "pocket",
                _build_method_panel(
                    "pocket",
                    "pocket",
                    lambda: gr.File(
                        label="Target PDB File",
                        file_types=[".pdb"],
                        type="filepath",
                    ),
                    visible=False,
                ),
            ),
            (
                "paired_pocket",
                _build_method_panel(
                    "paired_pocket",
                    "paired_pocket",
                    lambda: gr.File(
                        label="Target PDB File",
                        file_types=[".pdb"],
                        type="filepath",
                    ),
                    visible=False,
                ),
            ),
        ]

        switch_outputs = [
            *[panel for _name, panel in panels],
            *[button for _name, button in nav_buttons],
        ]
        for name, button in nav_buttons:
            button.click(
                fn=lambda selected=name: _switch_method(
                    selected,
                    panels,
                    nav_buttons,
                ),
                inputs=[],
                outputs=switch_outputs,
                api_name=f"show_{name}",
                show_progress="hidden",
            )

    return demo


demo = create_demo()


if __name__ == "__main__":
    demo.queue().launch(
        server_name=os.getenv("GRADIO_SERVER_NAME", "10.4.5.115"),
        server_port=int(os.getenv("GRADIO_SERVER_PORT", "7860")),
        show_error=True,
    )
