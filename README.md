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
* *Describe any prerequisites, libraries, OS version, etc., needed before installing program. ex. Windows 10*

We used [Python]( https://www.python.org/) as the main programming language and [GoogleColab](https://colab.research.google.com/) as the platform for code collaboration.

The main libraries and frameworks used in this project include:
* [DeepChem](https://deepchem.io/)
* [DGL](https://www.dgl.ai/)
* [Streamlit](https://streamlit.io/)
* [GitHub](https://github.com/)

### Installing

* How/where to download your program
* Any modifications needed to be made to files/folders
* 
All code work completed for the core components 1 & 2 of this project has been cleaned and saved as '.ipynb' files in the 'notebooks' folder in this [GitHub repository](https://github.com/luyingamypei/capstone_ui). To reproduce these parts of the project, one must install GoogleColab or another appropriate notebook environment and modify the relative directories to absolute paths work in their own environment before executing the code. To install a package within the notebook environment, use the code below as an example: 
```
!pip install deepchem
```

Our work on last core component has been integrated and stored in the same repository on Github, which is mandatory for publishing the user interface via Streamlit. This includes a main 'app.py' file that can deploy the app on the user interface via Streamlit Cloud along with other files such as required libraries and packages, pre-trained models and pre-built algorithms. Detailed instructions for execution are provided in the next section.

### Executing program

* *How to run the program
* *Step-by-step bullets
```
code blocks for commands
```

To reproduce our work in this project, please follow the steps below:
1. Exploratory Data Analysis (EDA): 
   For generating results and figures in EDA, execute the code in the 'Capstone_EDA.ipynb' file .
2. For training the toxicity prediction models and performing preliminary evaluations, run the code in 'Capstone_prelim_model_testing.ipynb' the toxicity model including model training and preliminary evaluations: 
4. Execute the code in 

## Help

Any advise for common problems or issues.
```
command to run if program contains helper info
```

## Authors

Contributors names and contact info

ex. Dominique Pizzie  
ex. [@DomPizzie](https://twitter.com/dompizzie)
Jonathan Gable @[gablejo@umich]()
Amy Pei @[Email] (luyingp@umich)

## Version History

* 0.2
    * Various bug fixes and optimizations
    * See [commit change]() or See [release history]()
* 0.1
    * Initial Release

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgments

Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)
