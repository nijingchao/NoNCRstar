For any questions about the source code and datasets, please contact [Jingchao Ni](http://nijingchao.github.io) (jingchao.ni@case.edu), thanks.

# Publication

* [Jingchao Ni](http://nijingchao.github.io), Mehmet Koyuturk, Hanghang Tong, Jonathan Haines, Rong Xu and Xiang Zhang. "Disease gene prioritization by integrating tissue-specific molecular networks using a robust multi-network model". BMC Bioinformatics, 2017.

# Requirements

* All codes are tested using MATLAB R2013a.


# Datasets

* P_G_NoSN_PPICenter.mat: the NoSN constructed using MAS threshold 0.4, disease similarity 0.4, tissue-specific gene expression level threshold 9. The tissue-specific PPI networks are the center networks.
* P_G_NoSN_GCNCenter.mat: the NoSN constructed using MAS threshold 0.4, disease similarity 0.4, tissue-specific gene expression level threshold 9. The center networks can be either tissue-specific gene co-expression networks or tissue-specific PPI networks.
* P_G_NoSN_PPICenter_2GCNs.mat: the NoSN constructed using MAS threshold 0.4, disease similarity 0.4, tissue-specific gene expression level threshold 9 and 300 for the first and second set of tissue-specific gene co-expression networks, respectively. The tissue-specific PPI networks are the center networks.

### Details

* PhenotypeSimNet: the disease similarity network.
* PhenotypeID: the MIM IDs of diseases in the disease similarity network.
* TSGeneNets: the tissue-specific molecular networks corresponding to the diseases in the disease similarity network. Each disease corresponds to at least one tissue-specific molecular network.
* TSGeneNetsID: the corresponding Entrez IDs of genes in TSGeneNets.
* AllGeneID: the Entrez IDs of all genes in all tissue-specific molecular networks.
* Seeds: the known causal genes of diseases in each tissue-specific molecular network, corresponding to TSGeneNets.
* TissueDict: the dictionary of the names of the most associated tissues of the diseases in the disease similarity network.

# Algorithm CR

* CR_CrossValidation.m: leave-one-out cross validation of CR.
* CR_Precomputation.m: the precomputation step of CR.
* CR: CR power method.
* J_CR.m: the objective function value of CrossRank.
* AUCEvaluation.m and AUCValue.m: AUC value evaluation with up to 50, 100, 300, 500, 700 and 1000 false positives.

Run CR_CrossValidation.m to see the evaluation results of the leave-one-out cross validation of CR algorithm. CR only uses the center networks in the NoSN datasets.

# Algorithm CRstar

* CRstar_CrossValidation.m: leave-one-out cross validation of CRstar.
* CRstar_Precomputation.m: the precomputation step of CRstar.
* CRstar: CRstar power method.
* J_CRstar.m: the objective function value of CrossRankStar.
* AUCEvaluation.m and AUCValue.m: AUC value evaluation with up to 50, 100, 300, 500, 700 and 1000 false positives.

Run CRstar_CrossValidation.m to see the evaluation results of the leave-one-out cross validation of CRstar algorithm.

# Algorithm WCRstar

* WCRstar_CrossValidation.m: leave-one-out cross validation of WCRstar.
* WCRstar_Precomputation.m: the precomputation step of WCRstar.
* WCRstar: WCRstar alternating minimization approach.
* J_WCRstar.m: the objective function value of Weighted CrossRankStar.
* J_aux.m: the center-auxiliary network inconsistency value of Weighted CrossRankStar.
* AUCEvaluation.m and AUCValue.m: AUC value evaluation with up to 50, 100, 300, 500, 700 and 1000 false positives.

Run WCRstar_CrossValidation.m to see the evaluation results of the leave-one-out cross validation of WCRstar algorithm. WCRstar uses P_G_NoSN_PPICenter_2GCNs.mat dataset as default.

# Note

* CR, CRstar and WCRstar have precomputation steps, which may take some time. The precomputation steps only need to be computed once for a dataset. If a precomputation file exists for a dataset, CR, CRstar and WCRstar will detect it and start ranking directly.