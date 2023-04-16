# Team Detox Capstone Project 

This project aims to build a web application for drug discovery scientists to input a molecule of interest, predict its toxicity and suggest new ideas for back-up candidates. 

## Description
*An in-depth paragraph about your project and overview of use.*

This project contains 3 core components:
1. Building a classification model for predicting the toxicity of a chemical molecule
2. Developing an idea generation tool for suggesting alternative molecules
3. Deploying the product (1 & 2) through an interactive user interface 

## Getting Started

### Dependencies
*Describe any prerequisites, libraries, OS version, etc., needed before installing program. ex. Windows 10*

We used [Python]( https://www.python.org/) as the main programming language and [GoogleColab](https://colab.research.google.com/) as the platform for code collaboration.

The main libraries and frameworks used in this project include:
* [DeepChem](https://deepchem.io/)
* [DGL](https://www.dgl.ai/)
* [Streamlit](https://streamlit.io/)
* [GitHub](https://github.com/)

### Installing

*How/where to download your program & any modifications needed to be made to files/folders*

All code work completed for the core components 1 & 2 of this project has been cleaned and saved as '.ipynb' files in the 'notebooks' folder in this [GitHub repository](https://github.com/luyingamypei/capstone_ui). To reproduce these parts of the project, one must install GoogleColab or another appropriate notebook environment and modify the relative directories to absolute paths work in their own environment before executing the code. If GoogleColab is used, Google Drive can be mounted to the runtime using the following code:
```
from google.colab import drive
drive.mount('/content/gdrive')
```
These notebooks contain code for installing all required libraries and packages for execution. To install any package within the notebook environment, use the code below as an example for installing the DeepChem library: 
```
!pip install deepchem
```

Our work on last core component has been integrated and stored in the same repository on Github, which is mandatory for publishing the user interface via Streamlit. This includes a main 'app.py' file that can deploy the app on the user interface via Streamlit Cloud along with other files such as required libraries and packages, pre-trained models and pre-built algorithms. To upgrade current packages or install additional packages for the app, modify or add them in the 'requirement.txt' file. Additionally, to enable any development work for the web application via Streamlit, one must sign up for an account on [Streamlit Cloud](https://streamlit.io/cloud) and set up a GitHub repository.


### Executing program

To reproduce our work in this project, please follow the steps below:
1. **Exploratory Data Analysis (EDA)**: 
   - For generating results and figures in EDA, execute all code cells in the 'Capstone_EDA.ipynb' file .
2. **Predictive Model Development**:
   - For training the toxicity prediction models and performing preliminary evaluations, run all code cells in 'Capstone_prelim_model_testing.ipynb'.
   - For performing the advanced model evaluations, execute all code cells in 'Capstone_advanced_model_evaluation.ipynb'.
3. **Idea Generation**:
   - For developing similarity search algorithm, refer to this notebook?
   - For developing idea tool via the matched molecular pairs approach, refer to this other notebook?  
4. **User interface**:
   - To test the experience on our fully deployed web-application: visit [https://luyingamypei-capstone-ui-app-amgnts.streamlit.app/](https://luyingamypei-capstone-ui-app-amgnts.streamlit.app/), 
   - To enable the development of this web application:
   
      a. Clone the GitHub repository with all required folders and files that can be found in the current [GitHub repository](https://github.com/luyingamypei/capstone_ui).  
      - app.py: main code for building the Streamlit application and deploying the product  
      - requirements.txt: contains all libraries and packages required for app.py to work
      - packages.txt: contains an extension package (libxrender1) required
      - 'AFPmodel' folder containing the pre-trained model stored as checkpoints 
      - mmpdb_ideas.py: pre-built algorithm for idea generation that will be loaded by app.py
      - AllHepG2.mmpdb: library of new molecules previously generated using the HepG2 dataset

      b. Log in to Streamlit Cloud and follow the instructions to deploy the app. For more information, refer to [How to Deploy an App on Streamlit Cloud](https://docs.streamlit.io/streamlit-community-cloud/get-started/deploy-an-app)

## Help

Common issues that may occur during program execution:
1. When cloning the GitHub repository, the 'HepG2.mmpdb'(395 MB) is a large file that may experience issues because GitHub has set an upper limit of 100 MB on the files one can upload. To resolve this issue, refer to [Git Large File Storage](https://git-lfs.com/), an open source Git extension for versioning large files.
2. When loading the Detox web application, a specific [ValueError](https://discuss.streamlit.io/t/valueerror-setting-an-array-element-with-a-sequence/40272) may occur randomly but would usually fix itself within hours, which seems to be a recent internal problem that should be addressed by the development team at Streamlit.     

## Authors

Jonathan Gable @Email: gablejo@umich.edu

Amy Pei @Email: luyingp@umich.edu

## Version History
* 0.2
    * Various bug fixes and optimizations
    * See [commit change]() or See [release history]()
* 0.1
    * Initial Release

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgments

*Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)
