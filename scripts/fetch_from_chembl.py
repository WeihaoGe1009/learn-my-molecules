from chembl_webresource_client.new_client import new_client
import requests
import yaml
from pathlib import Path

# Functional group → compound name mapping
group_to_compound = {
    "carboxylic acid": "acetic acid",
    "phenol": "phenol",
    "amine": "methylamine",
    "amide": "acetamide",
    "thiol": "methanethiol",
    "phosphate": "phosphoric acid",
    "aldehyde": "formaldehyde",
    "ketone": "acetone",
    "alcohol": "ethanol",
    "ether": "diethyl ether",
    "ester": "ethyl acetate",
    "nitrile": "acetonitrile",
    "guanidinium": "guanidine",
    "urea": "urea",
    "imine": "phenylhydrazone",
    "disulfide": "lipoic acid"
}

group_to_compound.update({
    "hydroxyl": "ethanol",                      # simple alcohol
    "carbonyl": "acetone",                      # ketone as representative carbonyl
    "enol": "acetylacetone",                    # tautomer of acetaldehyde
    "quaternary ammonium": "tetramethylammonium",  # stable quaternary amine
    "isoprenyl": "prenol",                    # base unit of isoprenoid
    "sulfonate": "methanesulfonic acid",        # common sulfonate group
    "azide": "CHEMBL3236174",                    # small organic azide
    "alkyne": "propyne",                        # terminal alkyne
    "boronic acid": "phenylboronic acid",       # common boronic acid in synthesis
    "maleimide": "maleimide",                   # standard compound, used directly
    "fluoro": "fluoromethane",                  # simplest fluoroalkane
    "nitro": "nitromethane",                    # simplest nitroalkane
})


# Load list of groups
input_file = Path("../data/functional_group_list.txt")

with input_file.open() as f:
    groups = [line.strip() for line in f if line.strip()]

output = {}

for group in groups:
    compound_name = group_to_compound.get(group.lower())
    if not compound_name:
        print(f"No mapping for: {group}")
        continue

    print(f"Processing: {group} → {compound_name}")
    results = {}

    # Search in ChEMBL
    if compound_name.upper().startswith("CHEMBL"):
        try:
            mol = new_client.molecule.get(compound_name)
            if mol:
                mols = [mol]
        except Exception as e:
            print(f"Failed to retrieve molecule for {compound_name}: {e}")
    else:
        mols = new_client.molecule.filter(pref_name__icontains=compound_name)
        if not mols:
            mols = new_client.molecule.filter(molecule_synonyms__molecule_synonym__iexact=compound_name)

    if not mols:
        print(f"No ChEMBL entry for: {compound_name}")
        continue

    mol = mols[0]
    chembl_id = mol["molecule_chembl_id"]
    results["representative_compound"] = compound_name
    results["chembl_id"] = chembl_id

    # SMILES
    structure = mol.get("molecule_structures", {})
    smiles = None
    for k, v in structure.items():
        if "smiles" in k.lower() and v:
            smiles = v
            break
    results["smiles"] = smiles


    # pKa fields
    pka_fields = {}
    for k, v in mol.get("molecule_properties", {}).items():
        if "pka" in k.lower() and v not in (None, "", "N/A"):
            try:
                pka_fields[k] = float(v)
            except ValueError:
                pka_fields[k] = v
    results["pKa"] = pka_fields or None

    # Target interactions
    bioactivities = new_client.activity.filter(molecule_chembl_id=chembl_id).only(
        "target_pref_name"
    )
    targets = {
        act["target_pref_name"]
        for act in bioactivities
        if act.get("target_pref_name")
           and not any(kw in act["target_pref_name"].lower() for kw in ["mus musculus", "no relevant", "homo sapiens"])
           and not act["target_pref_name"].strip().isdigit()
    }
    results["target_interactions"] = sorted(targets) if targets else None

    output[group] = results

# Save to YAML
with open("../data/functional_group_data.yaml", "w") as f:
    yaml.dump(output, f, sort_keys=False)

print("\n Done. Saved to ../data/functional_group_data.yaml")

