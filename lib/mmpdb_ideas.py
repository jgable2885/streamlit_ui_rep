import subprocess
import pandas as pd
from distutils.spawn import find_executable
from rdkit.Chem import PandasTools
from rdkit.Chem.Descriptors import qed
from rdkit import Chem
from rdkit import DataStructs


def generate_ideas(input_smi, prop="pCC50", database="AllHepG2.mmpdb", 
					delta="-1", max_ideas=None, output="out.tmp",
					sorting_weights=(0.7,0.3)):
    if find_executable("mmpdb") is None:
    	raise Exception("mmpdb is not installed")
    subprocess.run(["mmpdb", "transform", "-s", input_smi, 
                    database, "-p", prop, "--where", "avg<"+str(delta), "-o", output])
    ideas_df = pd.read_csv(output,sep="\t")
    PandasTools.AddMoleculeColumnToFrame(ideas_df,'SMILES','Molecule')
    ideas_df['qed'] = ideas_df['Molecule'].apply(qed)
    
    fps = [Chem.RDKFingerprint(mol) for mol in ideas_df['Molecule']]
    q_fp = Chem.RDKFingerprint(Chem.MolFromSmiles(input_smi))
    sims = DataStructs.BulkTanimotoSimilarity(q_fp, fps)
    ideas_df['sims'] = sims
    ideas_df['sort_by'] = sorting_weights[0]*ideas_df['sims']+sorting_weights[1]*ideas_df['qed']
    
    ideas_df.sort_values('sort_by', ascending=False, inplace=True)
    if max_ideas is not None:
    	return ideas_df[0:max_ideas]
    else:
    	return ideas_df
    
    