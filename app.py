import streamlit as st
import deepchem as dc
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from itertools import islice
from PIL import Image
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem import MolFromSmarts
from mmpdb_ideas import generate_ideas
import numpy as np
import altair as alt

def get_tox_class(prob_true):
	if prob_true > 0.6:
		return "Likely toxic"
	elif prob_true < 0.4:
		return "Likely non-toxic"
	else:
		return "Inconclusive"

def get_highlights(pattern, molecule):
    ##adapted from https://www.rdkit.org/docs/GettingStartedInPython.html
    if isinstance(molecule, str):
        mol = MolFromSmiles(molecule)
    elif isinstance(molecule, Chem.rdchem.Mol):
        mol = molecule
    else:
        raise TypeError('molecule must be either a smiles string or rdkit.Chem.rdchem.Mol object')
    if isinstance(pattern, str):
        patt = MolFromSmarts(pattern)
    else:
        raise TypeError('pattern must be a SMARTS or SMILES string')
    hit_ats = list(mol.GetSubstructMatch(patt))
    hit_bonds = []
    for bond in patt.GetBonds():
        aid1 = hit_ats[bond.GetBeginAtomIdx()]
        aid2 = hit_ats[bond.GetEndAtomIdx()]
        hit_bonds.append(mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())
    return hit_ats, hit_bonds

def find_tox21_match(query_smi, data='./data/tox21_y_fromDeepChem.csv', smiles_col='SMILES'):
	tox21_y_df = pd.read_csv(data)
	canonical_query = Chem.CanonSmiles(query_smi)
	if len(tox21_y_df[tox21_y_df[smiles_col]==canonical_query])>0:
		return True, tox21_y_df[tox21_y_df[smiles_col]==canonical_query]
	else:
		return False, []



st.write('Hello, welcome to the Detox App by Jonathan and Amy!')

input_smile = st.text_input('Please enter your compound of interest in SMILES format', 'SMILES Input')
	 
st.write('Press submit to have your molecule displayed below')
button1 = st.button('Submit button 1')

if st.session_state.get('button') != True:
	st.session_state['button'] = button1 # Saved the state

if st.session_state['button'] == True:
	st.write("button1 is True")
	st.write("Your compound of interest is: ", input_smile)
	
	mol = Chem.MolFromSmiles(input_smile)

	filename = "%s%d.png" % ("test", 0)

	Draw.MolToFile(mol, filename)
	st.image(Image.open(filename),caption='Compound structure')


	tox21_tasks3, tox21_datasets3, transformers3 = dc.molnet.load_tox21(featurizer=dc.feat.MolGraphConvFeaturizer(use_edges=True))
	st.write("Featurized!")
	model_reload=dc.models.AttentiveFPModel(n_tasks=12, batch_size=50, mode='classification',learning_rate=0.001, random_state=2, model_dir='AFPmodel')
	model_reload.restore()
	st.write("Model reloaded!")

	smiles = [input_smile]
	featurizer3 = dc.feat.MolGraphConvFeaturizer(use_edges=True)
	new_smile3 = featurizer3.featurize(smiles)
	
	## Check whether input smiles exactly matches molecule in Tox21 dataset
	match_found, data = find_tox21_match(smiles[0])
	if match_found == True:
		st.write('Exact match found, returning experimental results not predictions')
		exp_df = data.T.drop('SMILES')
		exp_df.reset_index(inplace=True)
		exp_df.columns = ['Assay Name', 'Prob Tox']
		exp_df['Probability of Toxicity'] = exp_df['Prob Tox'].astype(float).map(lambda n: '{:.2%}'.format(n))
		exp_df['tox_class'] = exp_df['Prob Tox'].apply(get_tox_class)
		exp_df.rename(columns={'tox_class':'Toxicity Class'}, inplace=True)
		
		st.table(exp_df.loc[:,['Assay Name', 'Toxicity Class', 'Probability of Toxicity']])
		chart_data = pd.DataFrame(exp_df[['Assay Name','Toxicity Class','Prob Tox','Probability of Toxicity']], 
								  columns=['Assay Name','Toxicity Class', 'Prob Tox','Probability of Toxicity'])

	else:
	
		preds = model_reload.predict_on_batch(new_smile3, transformers=transformers3)
		preds_df = pd.DataFrame(preds[0], columns=['Prob False','Prob Tox'])
		preds_df['Assay Name'] = tox21_tasks3
		preds_df['Probability of Toxicity'] = preds_df['Prob Tox'].astype(float).map(lambda n: '{:.2%}'.format(n))
		preds_df['tox_class'] = preds_df['Prob Tox'].apply(get_tox_class)
		preds_df.rename(columns= {'tox_class':'Toxicity Class'}, inplace = True)
	
		st.table(preds_df.loc[:,['Assay Name','Toxicity Class', 'Probability of Toxicity']])
		input_tox_count = np.sum(preds_df['Prob Tox'] > 0.6)
		st.write("{} predictions suggest toxicity".format(input_tox_count))
		chart_data = pd.DataFrame(preds_df[['Assay Name','Toxicity Class','Prob Tox','Probability of Toxicity']], 
								  columns=['Assay Name','Toxicity Class', 'Prob Tox','Probability of Toxicity'])

	chart = alt.Chart(chart_data, title="Tox Predictions for Your Input Molecule").mark_bar(color='darkred').encode(
		x = alt.X('Prob Tox:Q', scale=alt.Scale(domain=[0, 1.0]), axis=alt.Axis(format='%'), title='Probability of Toxicity'),
		y='Assay Name:N', 
		tooltip=['Assay Name','Toxicity Class', 'Probability of Toxicity'])
	st.altair_chart(chart, use_container_width=True) 
	
	ideas_df = generate_ideas(input_smile, database="AllHepG2.mmpdb")
	ideas_df.reset_index(inplace=True)

	featurized_ideas = featurizer3.featurize(ideas_df['SMILES'])
	idea_preds = model_reload.predict_on_batch(featurized_ideas, transformers=transformers3)
	toxicity_counts = [np.sum(pred>0.6) for pred in idea_preds[:,:,1]]
	ideas_df['toxicity_counts'] = toxicity_counts
	ideas_df.sort_values(by=['toxicity_counts','sort_by'], ascending=[True,False], inplace=True)
	genre = st.radio("How many new ideas would you like to generate?",('None', '1', '3', '5'))
	
	if st.button('Submit button 2'):
		st.write("Do your logic here")
		if genre == 'None':
			st.write("You have opted out for alternative candidate generation.")
		elif genre == '5':
			st.write("Let's see ", genre, " alternative ideas!")		
			idea_mols = [Chem.MolFromSmiles(smi) for smi in ideas_df['SMILES'][0:5]]
			idea_patts = [p for p in ideas_df['pCC50_to_smiles']]
			idea_legends = ["Predicted toxicities: {}".format(count) for count in ideas_df['toxicity_counts']]
			idea_hit_atoms = []
			idea_hit_bonds = []
			for p,m in zip(idea_patts, idea_mols):
				hit_ats, hit_bonds = get_highlights(p,m)
				idea_hit_atoms.append(hit_ats)
				idea_hit_bonds.append(hit_bonds)
			grid_img = Draw.MolsToGridImage(idea_mols, legends=idea_legends[0:5], highlightAtomLists=idea_hit_atoms, highlightBondLists=idea_hit_bonds)
			grid_img.save("tmp_grid.png")
			st.image(Image.open("tmp_grid.png"), caption='Idea Structures')
		elif genre == '3':
			st.write("Let's see ", genre, " alternative ideas!")
			idea_mols = [Chem.MolFromSmiles(smi) for smi in ideas_df['SMILES'][0:3]]
			idea_patts = [p for p in ideas_df['pCC50_to_smiles']]
			idea_legends = ["Predicted toxicities: {}".format(count) for count in ideas_df['toxicity_counts']]
			idea_hit_atoms = []
			idea_hit_bonds = []
			for p,m in zip(idea_patts, idea_mols):
				hit_ats, hit_bonds = get_highlights(p,m)
				idea_hit_atoms.append(hit_ats)
				idea_hit_bonds.append(hit_bonds)
			grid_img = Draw.MolsToGridImage(idea_mols, legends=idea_legends[0:3], highlightAtomLists=idea_hit_atoms, highlightBondLists=idea_hit_bonds)
			grid_img.save("tmp_grid.png")
			st.image(Image.open("tmp_grid.png"), caption='Idea Structures')
		elif genre == '1':
			st.write("Let's see ", genre, " alternative ideas!")
			idea_mols = [Chem.MolFromSmiles(smi) for smi in ideas_df['SMILES'][0:1]]
			idea_patts = [p for p in ideas_df['pCC50_to_smiles']]
			idea_legends = ["Predicted toxicities: {}".format(count) for count in ideas_df['toxicity_counts']]
			idea_hit_atoms = []
			idea_hit_bonds = []
			for p,m in zip(idea_patts, idea_mols):
				hit_ats, hit_bonds = get_highlights(p,m)
				idea_hit_atoms.append(hit_ats)
				idea_hit_bonds.append(hit_bonds)
			grid_img = Draw.MolsToGridImage(idea_mols, legends=idea_legends[0:1], highlightAtomLists=idea_hit_atoms, highlightBondLists=idea_hit_bonds)
			grid_img.save("tmp_grid.png")
			st.image(Image.open("tmp_grid.png"), caption='Idea Structures')
