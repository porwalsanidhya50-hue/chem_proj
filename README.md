
# ⚗️ Computational Chemistry Explorer

### 🧪 Explore Molecules • Visualize 3D Structures • Understand Chemical Properties

A Python-based interactive web application built with **Streamlit** that allows users to explore chemical compounds using their chemical formulas. It integrates **PubChem** for compound information and **RDKit** for molecular structure generation, visualization, and basic chemical analysis.

The project demonstrates how computational chemistry, cheminformatics, and data visualization can be combined in a simple, user-friendly application.

---

## ✨ Features

- 🔍 **Chemical Formula Search** — Enter a chemical formula to search for matching compounds.
- 🌐 **PubChem Integration** — Retrieve available compound information from PubChem.
- 🧬 **Molecular Structure Visualization** — Generate and display 2D molecular structures.
- 🧊 **3D Molecular Visualization** — Generate 3D conformers and explore molecular geometry interactively.
- ⚖️ **Molecular Weight Analysis** — Calculate or retrieve molecular mass.
- 🧩 **Elemental Composition** — Parse formulas and identify their constituent elements.
- 📊 **Chemical Properties** — Display available properties such as heavy atom count and rotatable bonds.
- ⚡ **Approximate Dipole Estimation** — Estimate a molecular dipole from calculated partial charges and 3D coordinates, where supported.
- 📋 **Compound Comparison** — Compare available properties of selected compounds.
- 🛡️ **Error Handling** — Handle invalid formulas, missing records, and unsuccessful lookups.

> ⚠️ The application is an educational prototype. Calculated properties are approximations where indicated and should not be treated as experimentally verified results.

---

## 🛠️ Technologies Used

| Technology | Purpose |
|---|---|
| 🐍 Python | Core programming language |
| 🎈 Streamlit | Interactive web interface |
| 🧬 RDKit | Molecular parsing, conformer generation, and chemical descriptors |
| 🔢 NumPy | Numerical calculations |
| 🐼 Pandas | Tabular data handling |
| 📈 Plotly | Interactive 3D visualization |
| 🌐 PubChem PUG REST API | Chemical compound data retrieval |

---

## 📁 Project Structure

```text
computational-chemistry-explorer/
│
├── app.py                 # Main Streamlit application
├── requirements.txt       # Python dependencies
├── README.md              # Project documentation
└── .gitignore             # Files excluded from Git
```

---

## 🚀 Getting Started

### 1️⃣ Prerequisites

Make sure you have installed:

- Python 3.10 or a compatible version
- pip
- An internet connection for PubChem searches

### 2️⃣ Clone the Repository

```bash
git clone <YOUR_REPOSITORY_URL>
cd computational-chemistry-explorer
```

If you have not uploaded the project to GitHub yet, simply open a terminal in your project folder.

### 3️⃣ Create a Virtual Environment (Recommended)

**Windows:**

```bash
python -m venv .venv
.venv\Scripts\activate
```

**macOS / Linux:**

```bash
python3 -m venv .venv
source .venv/bin/activate
```

### 4️⃣ Install Dependencies

```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
```

### 5️⃣ Run the Application

```bash
streamlit run app.py
```

The application should open in your default browser. If it does not, open the local URL displayed in your terminal.

---

## 📦 Dependencies

The `requirements.txt` file contains:

```text
streamlit
numpy
pandas
plotly
rdkit
```

---

## 🧑‍🔬 How to Use

1. Launch the application.
2. Enter a chemical formula, such as `H2O`, `CO2`, or `C6H6`.
3. Allow the application to search for matching compounds.
4. Select the intended compound from the available results.
5. Explore its 2D structure and interactive 3D representation.
6. Review its molecular weight and other available chemical properties.
7. Compare compounds where the application supports comparison.

### 🧪 Example Compounds

| Compound | Chemical Formula |
|---|---|
| 💧 Water | H₂O |
| 🌫️ Carbon dioxide | CO₂ |
| 🔥 Methane | CH₄ |
| 🧴 Ethanol | C₂H₆O |
| 💎 Benzene | C₆H₆ |
| 🧂 Potassium alum | KAl(SO₄)₂·12H₂O |

**Note:** A chemical formula does not always identify a unique compound. For example, C₂H₆O can represent ethanol or dimethyl ether. Always select the correct compound from the search results.

Hydrated ionic salts such as potassium alum may require specialized handling because a molecular structure generated from a SMILES string does not represent the complete ionic crystal.

---

## 🔬 Scientific Methodology

The application combines several computational steps:

1. **Formula parsing:** Validate elemental symbols and calculate the elemental composition.
2. **Compound identification:** Query PubChem for records matching the submitted formula.
3. **Structure conversion:** Convert a selected compound's SMILES representation into an RDKit molecular object.
4. **3D conformer generation:** Generate an approximate 3D geometry using RDKit's conformer-generation methods.
5. **Geometry optimization:** Attempt force-field optimization when supported by the molecule.
6. **Property analysis:** Retrieve or calculate selected molecular descriptors.
7. **Visualization:** Present molecular structures and numerical results through an interactive interface.

---

## ⚠️ Limitations

- 🔎 **Ambiguous formulas:** Different compounds can share the same molecular formula.
- 🌐 **Internet dependency:** PubChem searches require network access.
- 🧊 **3D geometry:** Generated conformers are computational estimates, not necessarily the experimentally observed geometry.
- ⚡ **Dipole moment:** Partial-charge-based estimates are approximate and depend on the molecular geometry and charge model.
- 🧂 **Ionic crystals:** Salts and hydrated compounds may not be represented correctly by ordinary single-molecule workflows.
- 🧪 **Scientific accuracy:** The application is intended for learning and demonstration, not for replacing quantum-chemical calculations or experimental measurements.

---

## 🎯 Project Objectives

- Understand the practical applications of computational chemistry.
- Explore molecular representations and chemical descriptors using RDKit.
- Integrate a public chemical database API into a Python application.
- Visualize molecular structures using interactive graphics.
- Build a user-friendly scientific application using Streamlit.
- Demonstrate the relationship between programming, chemistry, and scientific data analysis.

---

## 🔮 Future Enhancements

- 🧮 Integrate quantum chemistry software such as **Psi4** for more rigorous molecular calculations.
- 💎 Add crystal-structure visualization for ionic solids and hydrated salts.
- 📐 Display bond lengths, bond angles, and molecular geometry.
- ⚡ Improve dipole-moment calculations and compare them with reliable reference values.
- 📚 Add chemical reaction visualization and molecular property prediction.
- 📥 Export molecular structures, descriptors, and analysis reports.
- 🧠 Explore machine-learning models for predicting selected chemical properties.

---

## 🎓 Academic Context

**Project Type:** Internal Assessment / Mini Project  
**Domain:** Computational Chemistry, Cheminformatics, and Scientific Computing  
**Language:** Python

This project demonstrates the use of computational tools to retrieve, represent, visualize, and analyze chemical compounds through an interactive web application.

---

## 📚 References

- [PubChem](https://pubchem.ncbi.nlm.nih.gov/) — Chemical compound database
- [PubChem PUG REST API](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) — Programmatic access to chemical data
- [RDKit Documentation](https://www.rdkit.org/docs/) — Cheminformatics toolkit
- [Streamlit Documentation](https://docs.streamlit.io/) — Web application framework
- [Plotly Documentation](https://plotly.com/python/) — Interactive data visualization

---

## 👨‍💻 Author

**Sanidhya Porwal**

*Computational Chemistry Mini Project*

---

⭐ If you find this project useful, consider giving the repository a star!
