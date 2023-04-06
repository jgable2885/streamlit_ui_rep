import streamlit as st
import deepchem as dc
import pandas as pd
from rdkit import Chem
from itertools import islice
from PIL import Image
from rdkit import DataStructs
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
        mol = Chem.MolFromSmiles(molecule)
    elif isinstance(molecule, Chem.rdchem.Mol):
        mol = molecule
    else:
        raise TypeError('molecule must be either a smiles string or rdkit.Chem.rdchem.Mol object')
    if isinstance(pattern, str):
        patt = Chem.MolFromSmarts(pattern)
    else:
        raise TypeError('pattern must be a SMARTS or SMILES string')
    hit_ats = list(mol.GetSubstructMatch(patt))
    hit_bonds = []
    for bond in patt.GetBonds():
        aid1 = hit_ats[bond.GetBeginAtomIdx()]
        aid2 = hit_ats[bond.GetEndAtomIdx()]
        hit_bonds.append(mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())
    return hit_ats, hit_bonds

def find_tox21_match(query_smi, data, smiles_col='SMILES'):
	tox21_y_df = data
	canonical_query = Chem.CanonSmiles(query_smi)
	if len(tox21_y_df[tox21_y_df[smiles_col]==canonical_query])>0:
		return True, tox21_y_df[tox21_y_df[smiles_col]==canonical_query]
	else:
		return False, []

@st.cache_resource
def restore_model():
	model_reload=dc.models.AttentiveFPModel(n_tasks=12, batch_size=50, mode='classification',learning_rate=0.001, random_state=2, model_dir='AFPmodel')
	model_reload.restore()
	return model_reload

@st.cache_data
def load_tox21_data():
	tox21_tasks3, tox21_datasets3, transformers3 = dc.molnet.load_tox21(featurizer=dc.feat.MolGraphConvFeaturizer(use_edges=True))
	return tox21_tasks3, tox21_datasets3, transformers3
	
@st.cache_data
def load_tox21_y_csv():
	return pd.read_csv('./data/tox21_y_fromDeepChem.csv')

@st.cache_data
def generate_fingerprints(data, smiles_col='SMILES'):
	mols = [Chem.MolFromSmiles(x) for x in data[smiles_col]]
	fps = [Chem.RDKFingerprint(x) for x in mols]
	return fps
	
def display_idea_grid(ideas_df, num_ideas):
	idea_mols = [Chem.MolFromSmiles(smi) for smi in ideas_df['SMILES'][0:num_ideas]]
	idea_patts = [p for p in ideas_df['pCC50_to_smiles']]
	idea_legends = ["Predicted toxicities: {}".format(c) for c in ideas_df['toxicity_counts']]
	idea_hit_atoms = []
	idea_hit_bonds = []
	for p,m in zip(idea_patts, idea_mols):
		hit_ats, hit_bonds = get_highlights(p,m)
		idea_hit_atoms.append(hit_ats)
		idea_hit_bonds.append(hit_bonds)
	grid_img = Chem.Draw.MolsToGridImage(idea_mols, legends=idea_legends[0:num_ideas], 
										 highlightAtomLists=idea_hit_atoms, highlightBondLists=idea_hit_bonds)
	grid_img.save("tmp_grid.png")
	st.image(Image.open("tmp_grid.png"), caption='Idea Structures')

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

	Chem.Draw.MolToFile(mol, filename)
	st.image(Image.open(filename),caption='Compound structure')


	tox21_tasks3, tox21_datasets3, transformers3 = load_tox21_data()
	st.write("Featurized!")
	model_reload = restore_model()
	#model_reload=dc.models.AttentiveFPModel(n_tasks=12, batch_size=50, mode='classification',learning_rate=0.001, random_state=2, model_dir='AFPmodel')
	#model_reload.restore()
	st.write("Model reloaded!")

	smiles = [input_smile]
	featurizer3 = dc.feat.MolGraphConvFeaturizer(use_edges=True)
	new_smile3 = featurizer3.featurize(smiles)
	
	tox21_y_data = load_tox21_y_csv()
	fps = generate_fingerprints(tox21_y_data)
	input_fp = Chem.RDKFingerprint(Chem.MolFromSmiles(input_smile))
	
	##Find compounds similar to input in tox21 set
	sims = DataStructs.BulkTanimotoSimilarity(input_fp, fps)
	tox21_y_data['sim'] = sims
	tox21_y_data['toxicity_counts'] = tox21_y_data.loc[:,tox21_tasks3].sum(axis=1)
	tox21_y_data.sort_values(by=['sim','toxicity_counts'], ascending=[False, True], inplace=True)
	sim_mols = [Chem.MolFromSmiles(x) for x in tox21_y_data['SMILES'][0:3]]
	sim_legends = ["Experimental toxicities: {}".format(c) for c in tox21_y_data['toxicity_counts'][0:3]]
	sim_grid = Chem.Draw.MolsToGridImage(sim_mols, legends=sim_legends)
	sim_grid.save("sim_grid.png")
	st.image(Image.open("sim_grid.png"), caption='Top 3 most similar less toxic compounds from Tox21 database')
	
# 	grid_img = Chem.Draw.MolsToGridImage(idea_mols, legends=idea_legends[0:num_ideas], 
# 										 highlightAtomLists=idea_hit_atoms, highlightBondLists=idea_hit_bonds)
# 	grid_img.save("tmp_grid.png")
# 	st.image(Image.open("tmp_grid.png"), caption='Idea Structures')
	
	
	#ideas_df.sort_values(by=['toxicity_counts','sort_by'], ascending=[True,False], inplace=True)
	
	#toxicity_counts = [np.sum(pred>0.6) for pred in idea_preds[:,:,1]]
	
	match_found, data = find_tox21_match(smiles[0], tox21_y_data)
	if match_found == True:
		st.write('Exact match found, returning experimental results not predictions')
		exp_df = data.T.drop(['SMILES', 'toxicity_counts', 'sim'])
		exp_df.reset_index(inplace=True)
		exp_df.columns = ['Assay Name', 'Prob Tox']
		exp_df['Probability of Toxicity'] = exp_df['Prob Tox'].astype(float).map(lambda n: '{:.2%}'.format(n))
		exp_df['tox_class'] = exp_df['Prob Tox'].apply(get_tox_class)
		exp_df.rename(columns={'tox_class':'Toxicity Class'}, inplace=True)
		
		st.table(exp_df.loc[1:,['Assay Name', 'Toxicity Class', 'Probability of Toxicity']])
		chart_data = pd.DataFrame(exp_df[['Assay Name','Toxicity Class','Prob Tox','Probability of Toxicity']], 
								  columns=['Assay Name','Toxicity Class', 'Prob Tox','Probability of Toxicity'])
		chart_data = chart_data.loc[1:,:]

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
			display_idea_grid(ideas_df, int(genre))	

		elif genre == '3':
			st.write("Let's see ", genre, " alternative ideas!")
			display_idea_grid(ideas_df, int(genre))	

		elif genre == '1':
			st.write("Let's see ", genre, " alternative ideas!")
			display_idea_grid(ideas_df, int(genre))	
