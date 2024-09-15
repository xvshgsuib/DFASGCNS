DFASGCNS：A Prognostic Model for Ovarian Cancer Prediction Based on Dual Fusion Channels and Stacked Graph Convolu-tion
=====================
1.Introduction
-------------
In this paper, we propose a prognostic model for ovarian cancer prediction named Dual Fusion Channels and Stacked Graph Convolutional Neural Network (DFASGCNS). The DFASGCNS utilizes dual fusion channels to learn feature representations of different omics data and the associations between samples. Stacked graph convolutional network is used to comprehensively learn the deep and intricate correlation networks present in multi-omics data, enhancing the model's ability to represent multi-omics data. An attention mechanism is introduced to allocate different weights to important features of different omics data, optimizing the feature representation of multi-omics data. Experimental results demonstrate that compared to existing methods, the DFASGCNS model exhibits significant advantages in ovarian cancer prognosis prediction and survival analysis. Kaplan-Meier curve analysis results indicate significant differences in the survival subgroups predicted by the DFASGCNS model, contributing to a deeper understanding of the pathogenesis of ovarian cancer and providing more reliable auxiliary diagnostic information for the prognosis assessment of ovarian cancer patients.

2.Requirements
----------
Install Torch 1.10.0 and Python 3.6.11.

3.Data 
----------
Download Ovarian Cancer in TCGA Database: https://portal.gdc.cancer.gov/projects/TCGA-OV

4.Train
-----------
Run  ```models.python```  to start the training process.
