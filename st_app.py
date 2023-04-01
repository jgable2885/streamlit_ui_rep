import streamlit as st
import deepchem as dc
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from itertools import islice
#from IPython.display import Image, display
from PIL import Image
from rdkit.Chem.Draw import MolsToGridImage
from mmpdb_ideas import generate_ideas
#from stqdm import stqdm
import numpy as np
import altair as alt

def get_tox_class(prob_true):
	if prob_true > 0.6:
		return "Likely toxic"
	elif prob_true < 0.4:
		return "Likely non-toxic"
	else:
		return "Inconclusive"
	
st.write('Hello, welcome to the Detox App by Jonathan and Amy!')

if 'stage' not in st.session_state:
    st.session_state.stage = 0

def set_stage(stage):
    st.session_state.stage = stage

#Use session_state for nested submit forms
st.button('First Button', on_click=set_stage, args=(1,))

if st.session_state.stage > 0:
    # Some code
    st.button('Second Button', on_click=set_stage, args=(2,))
if st.session_state.stage > 1:
    # More code, etc
    st.button('Third Button', on_click=set_stage, args=(3,))
if st.session_state.stage > 2:
    st.write('The end')
st.button('Reset', on_click=set_stage, args=(0,))

form = st.form(key='my-form')
input_smile = form.text_input('Please enter your compound of interest in SMILES format', 'SMILES Input')

submit = form.form_submit_button('Submit', on_click=set_stage, args=(1,))
st.write('Press submit to have your molecule displayed below')

if submit:
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

	preds = model_reload.predict_on_batch(new_smile3, transformers=transformers3)
	preds_df = pd.DataFrame(preds[0], columns=['Prob False','Prob Tox'])
	preds_df['Assay'] = tox21_tasks3

	preds_df['tox_class'] = preds_df['Prob Tox'].apply(get_tox_class)
	st.table(preds_df.loc[:,['Assay','tox_class', 'Prob Tox']])
	input_tox_count = np.sum(preds_df['Prob Tox'] > 0.6)
	st.write("{} predictions suggest toxicity".format(input_tox_count))

	chart_data = pd.DataFrame(preds_df[['Assay','tox_class','Prob Tox']], columns=['Assay','tox_class', 'Prob Tox'])

	chart = alt.Chart(chart_data, title="Tox Predictions for Your Input Molecule").mark_bar(color='darkred').encode(
		x = alt.X('Prob Tox:Q', scale=alt.Scale(domain=[0, 1.0])),
		y='Assay:N', 
		tooltip=['Assay','tox_class', 'Prob Tox'])
	st.altair_chart(chart, use_container_width=True) 
	
	ideas_df = generate_ideas(input_smile, database="AllHepG2.mmpdb")
	ideas_df.reset_index(inplace=True)

	featurized_ideas = featurizer3.featurize(ideas_df['SMILES'])
	idea_preds = model_reload.predict_on_batch(featurized_ideas, transformers=transformers3)
	toxicity_counts = [np.sum(pred>0.6) for pred in idea_preds[:,:,1]]
	ideas_df['toxicity_counts'] = toxicity_counts
	ideas_df.sort_values(by=['toxicity_counts','sort_by'], ascending=[True,False], inplace=True)
	
	if st.session_state.stage > 1:
		form2 = st.form(key='2nd-form', on_click=set_stage, args=(2,))
		genre = form2.radio("How many new ideas would you like to generate?",('None', '1', '2', '3'))

		submit2 = form2.form_submit_button('Submit choice')

		if submit2:
			if genre == '3':
				st.write("Let's see ", genre, " alternative ideas!")

				idea_mols = [Chem.MolFromSmiles(smi) for smi in ideas_df['SMILES'][0:3]]
				idea_legends = ["Predicted toxicities: {}".format(count) for count in ideas_df['toxicity_counts']]
				grid_img = Draw.MolsToGridImage(idea_mols, legends=idea_legends[0:3])
				grid_img.save("tmp_grid.png")
				st.image(Image.open("tmp_grid.png"), caption='Idea Structures')
			elif genre == '2':
				st.write("Let's see ", genre, " alternative ideas!")

				idea_mols = [Chem.MolFromSmiles(smi) for smi in ideas_df['SMILES'][0:2]]
				idea_legends = ["Predicted toxicities: {}".format(count) for count in ideas_df['toxicity_counts']]
				grid_img = Draw.MolsToGridImage(idea_mols, legends=idea_legends[0:2])
				grid_img.save("tmp_grid.png")
				st.image(Image.open("tmp_grid.png"), caption='Idea Structures')
			elif genre == '1':
				st.write("Let's see ", genre, " alternative ideas!")

				idea_mols = [Chem.MolFromSmiles(smi) for smi in ideas_df['SMILES'][0:1]]
				idea_legends = ["Predicted toxicities: {}".format(count) for count in ideas_df['toxicity_counts']]
				grid_img = Draw.MolsToGridImage(idea_mols, legends=idea_legends[0:1])
				grid_img.save("tmp_grid.png")
				st.image(Image.open("tmp_grid.png"), caption='Idea Structures')
