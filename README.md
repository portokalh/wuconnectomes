# wuconnectomes
dipy based connectomes
and statistical analyses
to be expanded to include traits
Wenlin Wu & Alex Badea

1.	Streamline reconstruction is done using either of the two scripts:
	– CSA: CSA_revised_eg_sk_ww_binary_ALLanimal.py
	– CSD: CSD_det-reco.py
takes about 25 mins/subject; radius =20

2. Connectivity matrix generation

- Save_connec.pyn or Save_connec.ipynb
-output is saved : /Users/alex/code/Wenlin/preprocessing/results/

3.	combine connectivity csv to tensor
/Users/alex/code/Wenlin/preprocessing/ConcentrateMatrix2Tensor.py /Users/alex/code/Wenlin/preprocessing/ConcentrateMatrix2Tensor.ipynb

4.	Tensor PCA: relies on /Users/alex/code/Matlab/popNet_HOPCA rom Dr Zhang
/Users/alex/code/Matlab/popNet_HOPCA/ourscripts/TNPCAScripts.m

takes as input the tensor files from above
go to lines 73 or 76 and uncomments to generate TNPCA for select genotypes or ages
change nme of output files in line 105
        saves 'U – subject mode','V – network mode','D-lambda for principal components','percent_store-top 3 reflects x pwercent of variation ','obj-obj value function','connectivity-as tensor'

TNPCA_Clustering_plot to save figs
TNPCA identify subnet and copy results to xls sheet for LDA (fisher_training)

5. register bundles and do stats on bundles using:


Bundle based registration-for loop- Gen3&4-2way-WithFA.ipynb

Bundle based registration-for loop- Old&Young-2way-WithFA.ipynb	

