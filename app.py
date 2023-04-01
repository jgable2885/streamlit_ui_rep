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

st.write('Hello, this is a Streamlit test')

form = st.form(key='my-form')
input_smile = form.text_input('Please enter your compound of interest in SMILES format', 'SMILES Input')
submit = form.form_submit_button('Submit')

st.write('Press submit to have your molecule displayed below')

if submit:
    st.write("Your compound of interest is: ", input_smile)

    mol = Chem.MolFromSmiles(input_smile)

    filename = "%s%d.png" % ("test", 0)

    Draw.MolToFile(mol, filename)
    st.image(Image.open(filename),caption='Compound structure')

    st.write("Your predicted tox scores are:")
    tox21_tasks3, tox21_datasets3, transformers3 = dc.molnet.load_tox21(featurizer=dc.feat.MolGraphConvFeaturizer(use_edges=True))
#    tox21_tasks3, tox21_datasets3, transformers3 = dc.molnet.load_tox21(featurizer=dc.feat.ConvMolFeaturizer())
#     model_reload = dc.models.GraphConvModel(len(tox21_tasks3), mode="classification", model_dir="model")
#     model_reload.restore()
    model_reload=dc.models.AttentiveFPModel(n_tasks=12, batch_size=50, mode='classification', 
										learning_rate=0.001, random_state=2, model_dir='AFPmodel')
    model_reload.restore()

    smiles = [input_smile]
    featurizer3 = dc.feat.MolGraphConvFeaturizer(use_edges=True)
    new_smile3 = featurizer3.featurize(smiles)

    preds = model_reload.predict_on_batch(new_smile3, transformers=transformers3)
    preds_df = pd.DataFrame(preds[0], columns=['Prob False','Prob Tox'])
    preds_df['Assay'] = tox21_tasks3
    #preds_df.set_index('Assay')
    def get_tox_class(prob_true):
      if prob_true > 0.6:
        return "Likely toxic"
      elif prob_true < 0.4:
        return "Likely non-toxic"
      else:
        return "Inconclusive"
    preds_df['tox_class'] = preds_df['Prob Tox'].apply(get_tox_class)
    st.table(preds_df.loc[:,['Assay','tox_class', 'Prob Tox']])
    input_tox_count = np.sum(preds_df['Prob Tox'] > 0.6)
    st.write("{} predictions suggest toxicity".format(input_tox_count))
    #st.table(pd.DataFrame(preds[0], columns=['Prob Tox','Prob False']))

    ideas_df = generate_ideas(input_smile, database="AllHepG2.mmpdb")
    ideas_df.reset_index(inplace=True)

    featurized_ideas = featurizer3.featurize(ideas_df['SMILES'])
    idea_preds = model_reload.predict_on_batch(featurized_ideas, transformers=transformers3)
    toxicity_counts = [np.sum(pred>0.6) for pred in idea_preds[:,:,1]]
    ideas_df['toxicity_counts'] = toxicity_counts
    ideas_df.sort_values(by=['toxicity_counts','sort_by'], ascending=[True,False], inplace=True)

    idea_mols = [Chem.MolFromSmiles(smi) for smi in ideas_df['SMILES'][0:3]]
    idea_legends = ["Predicted toxicities: {}".format(count) for count in ideas_df['toxicity_counts']]
    grid_img = Draw.MolsToGridImage(idea_mols, legends=idea_legends[0:3])
    grid_img.save("tmp_grid.png")
    st.image(Image.open("tmp_grid.png"), caption='Idea Structures')
