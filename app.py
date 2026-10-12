
import re
import json
import urllib.parse
import urllib.request
import urllib.error
from collections import Counter

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, rdMolDescriptors


# =========================================================
# CONFIGURATION
# =========================================================

st.set_page_config(
    page_title="Molecular Polarity Lab",
    page_icon="⚗️",
    layout="wide"
)

API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
MAX_RESULTS = 20

st.title("⚗️ Molecular Polarity Lab")
st.caption(
    "Formula parser • Compound discovery • Molecular visualization"
)

st.info(
    "Enter a molecular formula to find matching compounds. "
    "Because different structures can share the same formula, "
    "select a specific compound before analysing its geometry."
)


# =========================================================
# FORMULA PARSER
# Supports elements, nested parentheses, square brackets,
# and hydrate separators such as CuSO4·5H2O.
# =========================================================

PERIODIC_TABLE = Chem.GetPeriodicTable()


def parse_formula(formula):
    formula = formula.strip()

    if not formula or len(formula) > 150:
        raise ValueError("Enter a formula between 1 and 150 characters.")

    # Normalize Unicode subscripts and common hydrate separators.
    subscript_map = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
    formula = formula.translate(subscript_map)
    formula = formula.replace("·", ".").replace("⋅", ".")

    if any(ch.isspace() for ch in formula):
        raise ValueError("Remove spaces from the chemical formula.")

    if re.search(r"[^A-Za-z0-9()[\].]", formula):
        raise ValueError(
            "Use element symbols, numbers, parentheses, brackets "
            "and hydrate separators only."
        )

    def parse_segment(text):
        tokens = re.findall(
            r"[A-Z][a-z]?|\d+|[()\[\]]",
            text
        )

        if "".join(tokens) != text:
            raise ValueError("Invalid chemical formula syntax.")

        pos = 0

        def parse_group(closing=None):
            nonlocal pos
            counts = Counter()
            found = False

            while pos < len(tokens):
                token = tokens[pos]

                if token in (")", "]"):
                    if closing != token:
                        raise ValueError("Mismatched brackets.")
                    pos += 1
                    return counts

                if token in ("(", "["):
                    found = True
                    expected = ")" if token == "(" else "]"
                    pos += 1
                    inner = parse_group(expected)

                    multiplier = 1
                    if pos < len(tokens) and tokens[pos].isdigit():
                        multiplier = int(tokens[pos])
                        pos += 1

                    if multiplier < 1:
                        raise ValueError("Subscripts must be positive.")

                    for element, number in inner.items():
                        counts[element] += number * multiplier

                elif re.fullmatch(r"[A-Z][a-z]?", token):
                    found = True
                    element = token

                    if PERIODIC_TABLE.GetAtomicNumber(element) == 0:
                        raise ValueError(f"Unknown element: {element}")

                    pos += 1
                    multiplier = 1

                    if pos < len(tokens) and tokens[pos].isdigit():
                        multiplier = int(tokens[pos])
                        pos += 1

                    if multiplier < 1:
                        raise ValueError("Subscripts must be positive.")

                    counts[element] += multiplier

                else:
                    raise ValueError("Unexpected number or token.")

            if closing:
                raise ValueError("A closing bracket is missing.")

            if not found:
                raise ValueError("Empty formula group.")

            return counts

        counts = parse_group()

        if pos != len(tokens):
            raise ValueError("Could not parse the complete formula.")

        return counts

    total = Counter()

    # A dot separates hydrate components.
    parts = re.split(r"[.]", formula)

    for part in parts:
        if not part:
            raise ValueError("Empty hydrate component.")

        # Optional coefficient at the beginning of a component.
        match = re.match(r"^(\d+)(?=[A-Z[(])", part)

        coefficient = 1
        if match:
            coefficient = int(match.group(1))
            part = part[match.end():]

        if coefficient < 1:
            raise ValueError("Component coefficients must be positive.")

        component = parse_segment(part)

        for element, number in component.items():
            total[element] += number * coefficient

    if not total:
        raise ValueError("No elements were found.")

    if sum(total.values()) > 10000:
        raise ValueError("Formula is too large to process.")

    return total


def hill_formula(counts):
    """Convert atom counts to conventional Hill notation."""
    elements = list(counts)

    if "C" in counts:
        order = ["C"]
        if "H" in counts:
            order.append("H")
        order.extend(sorted(e for e in elements if e not in ("C", "H")))
    else:
        order = sorted(elements)

    return "".join(
        element + (str(counts[element]) if counts[element] != 1 else "")
        for element in order
    )


def formula_mass(counts):
    return sum(
        PERIODIC_TABLE.GetAtomicWeight(
            PERIODIC_TABLE.GetAtomicNumber(element)
        ) * count
        for element, count in counts.items()
    )


# =========================================================
# PUBCHEM API
# =========================================================

def api_get(path, timeout=20):
    url = f"{API}/{path}"

    request = urllib.request.Request(
        url,
        headers={"User-Agent": "MolecularPolarityLab/1.0"}
    )

    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def search_pubchem(formula):
    encoded = urllib.parse.quote(formula, safe="")

    data = api_get(
        f"compound/fastformula/{encoded}/cids/JSON"
        f"?MaxRecords={MAX_RESULTS}"
    )

    return data.get("IdentifierList", {}).get("CID", [])


def fetch_compound_properties(cids):
    if not cids:
        return []

    cid_text = ",".join(str(cid) for cid in cids)

    fields = (
        "MolecularFormula,MolecularWeight,"
        "CanonicalSMILES,IUPACName"
    )

    data = api_get(
        f"compound/cid/{cid_text}/property/{fields}/JSON"
    )

    return data.get("PropertyTable", {}).get("Properties", [])


# =========================================================
# MOLECULAR STRUCTURE
# =========================================================

def make_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError("RDKit could not interpret this structure.")

    return mol


def make_3d(mol):
    """Create a 3D conformer; return a hydrogen-containing molecule."""
    mol3d = Chem.AddHs(Chem.Mol(mol))

    params = AllChem.ETKDGv3()
    params.randomSeed = 42

    status = AllChem.EmbedMolecule(mol3d, params)

    if status != 0:
        raise ValueError(
            "Could not generate a 3D conformer for this structure."
        )

    # MMFF is attempted first; UFF is a fallback.
    if AllChem.MMFFHasAllMoleculeParams(mol3d):
        AllChem.MMFFOptimizeMolecule(mol3d, maxIters=300)
    elif AllChem.UFFHasAllMoleculeParams(mol3d):
        AllChem.UFFOptimizeMolecule(mol3d, maxIters=300)

    return mol3d


def draw_3d(mol3d):
    conf = mol3d.GetConformer()
    positions = np.array([
        list(conf.GetAtomPosition(i))
        for i in range(mol3d.GetNumAtoms())
    ])

    symbols = [
        atom.GetSymbol()
        for atom in mol3d.GetAtoms()
    ]

    colors = [
        "#EF4444" if symbol == "O"
        else "#3B82F6" if symbol == "N"
        else "#22C55E" if symbol == "Cl"
        else "#F59E0B" if symbol == "S"
        else "#94A3B8" if symbol == "C"
        else "#E5E7EB"
        for symbol in symbols
    ]

    fig = go.Figure()

    # Draw bonds using actual conformer coordinates.
    for bond in mol3d.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        fig.add_trace(go.Scatter3d(
            x=[positions[i, 0], positions[j, 0]],
            y=[positions[i, 1], positions[j, 1]],
            z=[positions[i, 2], positions[j, 2]],
            mode="lines",
            line=dict(color="#94A3B8", width=6),
            hoverinfo="skip",
            showlegend=False
        ))

    fig.add_trace(go.Scatter3d(
        x=positions[:, 0],
        y=positions[:, 1],
        z=positions[:, 2],
        mode="markers+text",
        text=symbols,
        textposition="top center",
        marker=dict(
            size=7,
            color=colors,
            line=dict(width=1, color="#64748B")
        ),
        hovertemplate="%{text}<extra></extra>",
        showlegend=False
    ))

    fig.update_layout(
        height=470,
        margin=dict(l=0, r=0, t=15, b=0),
        scene=dict(
            xaxis_title="X (Å)",
            yaxis_title="Y (Å)",
            zaxis_title="Z (Å)",
            aspectmode="data"
        )
    )

    return fig


# =========================================================
# APPROXIMATE PARTIAL-CHARGE DIPOLE
# This is a force-field/charge-model estimate, not a
# quantum-chemical dipole moment.
# =========================================================

def approximate_dipole(mol3d):
    """
    Estimate dipole from RDKit Gasteiger partial charges.

    The vector is sum(q_i * r_i). Its magnitude is converted
    to Debye using 4.8032 D per e·Å.

    This is approximate and depends on the charge model,
    conformer, and coordinate origin for charged molecules.
    """

    mol = Chem.Mol(mol3d)

    try:
        Chem.rdPartialCharges.ComputeGasteigerCharges(mol)
    except Exception:
        return None, None

    conf = mol.GetConformer()
    dipole = np.zeros(3)

    for atom in mol.GetAtoms():
        charge_text = atom.GetProp("_GasteigerCharge")

        try:
            charge = float(charge_text)
        except (ValueError, TypeError):
            return None, None

        if not np.isfinite(charge):
            return None, None

        position = np.array(
            list(conf.GetAtomPosition(atom.GetIdx()))
        )

        dipole += charge * position

    magnitude = float(np.linalg.norm(dipole))
    debye = magnitude * 4.8032

    return dipole, debye


# =========================================================
# FORMULA INPUT
# =========================================================

with st.form("formula_form"):
    entered_formula = st.text_input(
        "Chemical formula",
        value="H2O",
        placeholder="Examples: H2O, C6H6, CH3COOH, Ca(OH)2",
        help=(
            "Use element symbols and numbers. Parentheses and "
            "nested groups are supported. Hydrate separators "
            "such as CuSO4·5H2O are supported for formula parsing."
        )
    )

    submitted = st.form_submit_button(
        "Analyse formula",
        type="primary",
        use_container_width=True
    )


if submitted:
    st.session_state.pop("compound_records", None)
    st.session_state.pop("selected_formula", None)
    st.session_state.pop("selected_counts", None)

    try:
        counts = parse_formula(entered_formula)
        normalized_formula = hill_formula(counts)

        st.session_state.selected_formula = normalized_formula
        st.session_state.selected_counts = dict(counts)

        with st.spinner("Searching PubChem for matching compounds..."):
            cids = search_pubchem(normalized_formula)

            if cids:
                records = fetch_compound_properties(cids)
            else:
                records = []

        st.session_state.compound_records = records
        st.session_state.search_finished = True

    except urllib.error.HTTPError as exc:
        st.error(
            f"PubChem returned HTTP error {exc.code}. "
            "Try again later or check your internet connection."
        )
    except urllib.error.URLError:
        st.error(
            "Could not connect to PubChem. Check your internet "
            "connection and try again."
        )
    except (ValueError, KeyError, json.JSONDecodeError) as exc:
        st.error(f"Could not process this formula: {exc}")


# =========================================================
# RESULTS
# =========================================================

if "selected_formula" in st.session_state:

    counts = st.session_state.selected_counts
    normalized_formula = st.session_state.selected_formula

    st.markdown("---")
    st.subheader("Formula Analysis")

    col1, col2, col3 = st.columns(3)

    col1.metric("Normalized Formula", normalized_formula)
    col2.metric("Number of Atoms", sum(counts.values()))
    col3.metric("Approx. Formula Mass", f"{formula_mass(counts):.3f} g/mol")

    composition = pd.DataFrame([
        {
            "Element": element,
            "Atom Count": count,
            "Atomic Mass (u)": round(
                PERIODIC_TABLE.GetAtomicWeight(
                    PERIODIC_TABLE.GetAtomicNumber(element)
                ), 4
            ),
            "Mass Contribution (g/mol)": round(
                PERIODIC_TABLE.GetAtomicWeight(
                    PERIODIC_TABLE.GetAtomicNumber(element)
                ) * count, 4
            )
        }
        for element, count in counts.items()
    ])

    st.write("Elemental composition")
    st.dataframe(composition, hide_index=True, use_container_width=True)

    records = st.session_state.get("compound_records", [])

    st.markdown("---")
    st.subheader("Matching Compounds")

    if not records:
        st.warning(
            "No matching compound was found in PubChem, or the "
            "search returned no usable records. The formula and "
            "composition are still available above. This does not "
            "prove that the substance does not exist."
        )

    else:
        st.caption(
            f"Showing up to {len(records)} records from PubChem. "
            "The full list may contain additional compounds."
        )

        options = {}

        for record in records:
            cid = record.get("CID")
            name = record.get("IUPACName") or "Name unavailable"
            smiles = record.get("ConnectivitySMILES") or record.get(
                "CanonicalSMILES"
            )

            if cid is None or not smiles:
                continue

            label = f"CID {cid} — {name}"
            options[label] = record

        if options:
            selected_label = st.selectbox(
                "Select the specific compound",
                list(options.keys())
            )

            selected = options[selected_label]
            smiles = (
                selected.get("ConnectivitySMILES")
                or selected.get("CanonicalSMILES")
            )

            try:
                mol = make_molecule(smiles)
                mol3d = make_3d(mol)

                st.markdown("---")
                st.subheader("Molecular Properties")

                c1, c2, c3 = st.columns(3)

                c1.metric(
                    "Molecular Weight",
                    f"{Descriptors.MolWt(mol):.3f} g/mol"
                )
                c2.metric(
                    "Heavy Atoms",
                    mol.GetNumHeavyAtoms()
                )
                c3.metric(
                    "Rotatable Bonds",
                    rdMolDescriptors.CalcNumRotatableBonds(mol)
                )

                left, right = st.columns([1, 1.4])

                with left:
                    st.subheader("2D Structure")
                    image_mol = Chem.RemoveHs(Chem.Mol(mol))
                    st.image(
                        Draw.MolToImage(
                            image_mol,
                            size=(450, 350)
                        ),
                        use_container_width=True
                    )

                with right:
                    st.subheader("3D Molecular Structure")
                    st.plotly_chart(
                        draw_3d(mol3d),
                        use_container_width=True
                    )

                st.markdown("---")
                st.subheader("Polarity Estimate")

                dipole_vector, dipole_debye = approximate_dipole(mol3d)

                if dipole_vector is None:
                    st.warning(
                        "Could not calculate a Gasteiger-charge "
                        "dipole estimate for this structure."
                    )
                else:
                    d1, d2 = st.columns(2)

                    d1.metric(
                        "Estimated Dipole Magnitude",
                        f"{dipole_debye:.3f} D"
                    )

                    d2.write("Estimated vector (Debye)")
                    d2.code(
                        "(" + ", ".join(
                            f"{x * 4.8032:.3f}"
                            for x in dipole_vector
                        ) + ")"
                    )

                    st.warning(
                        "This is a rough estimate from Gasteiger "
                        "partial charges and one generated conformer. "
                        "It is not an experimental or quantum-chemical "
                        "dipole moment. Do not use a universal numerical "
                        "threshold here to classify all molecules."
                    )

                st.markdown("### Compound Details")

                st.write("**PubChem CID:**", selected.get("CID"))
                st.write(
                    "**IUPAC name:**",
                    selected.get("IUPACName", "Not available")
                )
                st.write(
                    "**Canonical SMILES:**",
                    selected.get("ConnectivitySMILES")
                    or selected.get("CanonicalSMILES", "Not available")
                )

                cid = selected.get("CID")
                if cid:
                    st.markdown(
                        f"[Open compound record on PubChem]"
                        f"(https://pubchem.ncbi.nlm.nih.gov/compound/{cid})"
                    )

            except (ValueError, RuntimeError, KeyError) as exc:
                st.error(f"Could not analyse the selected structure: {exc}")

        else:
            st.warning(
                "PubChem returned records, but none contained a "
                "usable molecular structure."
            )


# =========================================================
# FOOTER
# =========================================================

with st.expander("Scientific notes and limitations"):
    st.markdown("""
    - A molecular formula specifies atom counts, not unique connectivity.
    - PubChem formula searches may return multiple structures.
    - Formula parsing and formula mass calculation work independently
      of whether a compound is found online.
    - RDKit creates and optimizes a molecular conformer using available
      force-field methods; this is not a quantum-chemical calculation.
    - Gasteiger charges provide an approximate partial-charge model.
      The resulting dipole estimate can differ from experimental values.
    - Salts, coordination compounds, unusual valence states, isotopes,
      and some charged structures may need specialized treatment.
    - Hydrate formula parsing calculates total elemental composition;
      it does not guarantee that PubChem treats a hydrate as one
      covalently connected molecule.
    """)

st.caption(
    "Molecular Polarity Lab • Python • Streamlit • RDKit • PubChem"
)
