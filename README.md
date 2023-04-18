# Team Detox Capstone Project 

A web application for predictive models and idea generation to improve preclinical toxicology. 

## Description
*An in-depth paragraph about your project and overview of use.*

This project aims to build a web application that would allow drug discovery scientists to input a molecule of interest, predict its toxicity and suggest new ideas for back-up candidates. 
There are 3 core components:
1. Building a classification model for predicting the toxicity of a chemical molecule
2. Developing an idea generation tool for suggesting alternative molecules
3. Deploying the product (1 & 2) through an interactive user interface 

## Getting Started

### Dependencies

We used [Python]( https://www.python.org/) as the main programming language and [GoogleColab](https://colab.research.google.com/) as the platform for code collaboration.

The main libraries and frameworks used in this project include:
* [DeepChem](https://deepchem.io/)
* [RDKit](https://www.rdkit.org/)
* [mmpdb](https://github.com/rdkit/mmpdb)
* [DGL](https://www.dgl.ai/) 
* [scikit-learn](https://scikit-learn.org/stable/)
* [Streamlit](https://streamlit.io/)
* [GitHub](https://github.com/)

### Installing

*How/where to download your program & any modifications needed to be made to files/folders*

All code work completed for the core components 1 & 2 of this project has been cleaned and saved as '.ipynb' files in the ['notebooks' folder](https://github.com/luyingamypei/capstone_ui/tree/main/notebooks) in this current GitHub repository. To reproduce these parts of the project, one must install GoogleColab or another appropriate notebook environment and either recreate the directory structure of this repo or modify the paths in the code to work in their own environment before executing the code. If GoogleColab is used, Google Drive can be mounted to the runtime using the following code:
```
from google.colab import drive
drive.mount('/content/gdrive')
```
The code for installing all the libraries and packages required for running the notebooks can be found within the notebooks. To install any package within the notebook environment, use the code below as an example for installing the DeepChem library: 
```
!pip install deepchem
```

Our work on last core component has been integrated and stored in the same repository on Github, which is mandatory for publishing the user interface via Streamlit. This includes a main 'app.py' file that can deploy the app on the user interface via Streamlit Cloud along with other files such as required libraries and packages, pre-trained models and pre-built algorithms. To upgrade current packages or install additional packages for the app, modify or add them in the 'requirement.txt' file. Additionally, to enable any development work for the web application via Streamlit, one must sign up for an account on [Streamlit Cloud](https://streamlit.io/cloud) and set up a [GitHub repository](https://docs.github.com/en/get-started/quickstart/create-a-repo).


### Executing program

To reproduce our work in this project, please follow the steps below:
1. **Exploratory Data Analysis (EDA)**: 
   - For generating results and figures in EDA, execute all code cells in the 'EDA_TrainingSetOnly.ipynb' file in the 'notebooks' folder.
   - For generating a CSV file containing Tox21 experimental data with canonical SMILES, refer to 'GenerateCSVwithTox21CanonicalSmilesAndExperimentalResults.ipynb' in the 'notebooks' folder.  
2. **Predictive Model Development**:
   - For training the toxicity prediction models and performing preliminary evaluations, run all code cells in 'ToxModelTrainingAndPreliminaryEvaluations.ipynb' in the 'notebooks' folder.
   - For performing the advanced model evaluations, execute all code cells in 'ToxModelAdvancedEvaluations.ipynb' in the 'notebooks' folder.
   - All of the trained deep learning models (except for the Weave model due to size restrictions) have been stored in the 'models' folder.
3. **Idea Generation**:
   - To retrieve the ChEMBL HepG2 data used for MMP database creation and some figure generation run all cells of the 'RetrieveChEMBLData.ipynb' file. Retrieving the ~58K rows from ChEMBL takes several hours. A CSV file, hepg2_output_df.csv, is available in the data directory of this repo.
   - For generation of MMP databases via the mmpdb package run all cells of the CreateMMPdatabases.ipynb file. This takes ~1 hour to execute and creates large files: ~400MB for the HepG2 MMP database and ~2.6GB for the Tox21 MMP database. The HepG2 database is contained in 'AllHepG2.mmpdb' in this repo.
   - A helper function that takes a SMILES string as input and returns a Pandas dataframe of ideas via the mmpdb package can be found in the mmpdb_ideas.py file in the lib directory.
   - Code for similarity search to retrieve the most similar Tox21 molecules to the user input can be found in 'app.py'
   - For evaluation of the idea generation approach and associated figures execute all cells of 'IdeaGenerationEvaluation.ipynb' in the 'notebooks' folder.  
4. **User interface**:
   - To test the experience on our fully deployed web-application: visit [https://luyingamypei-capstone-ui-app-amgnts.streamlit.app/](https://luyingamypei-capstone-ui-app-amgnts.streamlit.app/), 
   - To enable the development of this online application:
   
      a. Clone the GitHub repository with all required folders and files that can be found in the current [GitHub repository](https://github.com/luyingamypei/capstone_ui).  
      - app.py: main code for building the Streamlit application and deploying the product  
      - requirements.txt: contains all libraries and packages required for app.py to work
      - packages.txt: contains an extension package (libxrender1) required
      - 'AFPmodel' folder containing the pre-trained model stored as checkpoints that will be loaded by app.py
      - mmpdb_ideas.py: pre-built algorithm for idea generation that will be loaded by app.py
      - AllHepG2.mmpdb: library of new molecules previously generated using the HepG2 dataset

      b. Log in to Streamlit Cloud and follow the instructions to deploy the app. For more information, refer to [How to Deploy an App on Streamlit Cloud](https://docs.streamlit.io/streamlit-community-cloud/get-started/deploy-an-app)

## Help

Common issues that may occur during program execution:
1. When cloning the GitHub repository, the 'HepG2.mmpdb'(395 MB) is a large file that may experience issues because GitHub has set an upper limit of 25 MB per file on the files uploaded to a repository via a browser. To resolve this issue, refer to [Git Large File Storage](https://git-lfs.com/), an open source Git extension for versioning large files.
2. When loading the Detox web application, a specific [ValueError](https://discuss.streamlit.io/t/valueerror-setting-an-array-element-with-a-sequence/40272) may occur randomly but would usually fix itself within hours, which seems to be a recent internal problem that should be addressed by the development team at Streamlit.     

## Authors

Jonathan Gable: gablejo@umich.edu

Amy Pei: luyingp@umich.edu

Project Link: https://github.com/luyingamypei/capstone_ui

## Version History
    
* 0.1
    * Initial Release
    * See [commit change](https://github.com/luyingamypei/capstone_ui/commits/main) or See [release history](https://github.com/luyingamypei/capstone_ui/releases)

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details

## Acknowledgments

*Inspiration, code snippets, etc.*

* We would like to acknowledge the teams who developed the main libraries, tools and platforms that were crucial in our project, including but not limited to DeepChem, Streamlit, GitHub, Python, and GoogleColab.
* We would also like to thank everyone who contributed their inputs to our qualitative data inquiries which led to valuable insights for our data analysis.
* Finally, we would like to acklowdge and thank the entire MADS Capstone instructional team especially Dr. O'Brien, Michelle LeBlanc and all others who provided us with guidance and help throughout this project.


