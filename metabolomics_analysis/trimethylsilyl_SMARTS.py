import pandas as pd 
from rdkit import Chem
from rdkit.Chem import AllChem 

def revert_tms_derivatization(silylated):
    mol = Chem.MolFromSmiles(silylated)

    # Reaction SMARTS: [O,N,S:1]-Si(C)(C)C >> [O,N,S:1]-H
    # This identifies the heteroatom (mapped to :1) and strips the silyl group.
    rxn_smarts = '[O,N,S:1][Si](C)(C)C>>[*:1]'
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)

    # Run the reaction (this may produce multiple products if multiple TMS groups exist)
    products = rxn.RunReactants((mol,))
    
    if products:
        # Get the first result and convert to SMILES
        final_mol = products[0][0]
        Chem.SanitizeMol(final_mol) # Ensure chemical validity
        return Chem.MolToSmiles(final_mol)
    
    return silylated

