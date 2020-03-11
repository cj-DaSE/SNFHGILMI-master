
Type: MATLAB Code

Title: Heterogeneous graph inference based on similarity network fusion for predicting lncRNA-miRNA interaction

=================

Description: This program can be opened by matlab R2015b. After opening it, you can run it directly by running the script file named main.m.
This program implements the SNFHGILMI algorithm predicting lncRNA-miRNA 
interactions.

Files description:

1.Data: 
  1) Known_lncRNA_miRNA_association.mat: The id pairs of known lncRNA-miRNA interactions.

     lncRNA_id.txt: The id list of lncRNA name.

     miRNA_id.txt: The id list of miRNA name.

  2) Gaussian_kernel_data

      Gaussksim.m: Construct lncRNA/miRNA Gaussian interaction profile kernel similarity matrix.

      lnc_gaussian_similarity_matrix.mat: The lncRNA-lncRNA similarity matrix of Gaussian interaction profile.

      mi_gaussian_similarity_matrix.mat: The miRNA-miRNA similarity matrix of Gaussian interaction profile. 

   3) Sequence_data

     invalid_lnc_seq.mat: The ids of lncRNAs without lncRNA sequence information.

     lnc_seq_similarity_matrix.mat: The lncRNA-lncRNA similarity matrix of sequence.

     mi_seq_similarity_matrix.mat: The miRNA-miRNA similarity matrix of sequence.


2.Code:

  1) 5-foldcv

    main.m: The main program for the implemention of 5-fold cross validation.

    Data_fusion_SNF.m: Integrate the similarity matrix of lncRNA/miRNA.

    Random_roder_5fold.m: Generating the random sample lists.
   
    HGILMA.m: The subprogram for the implemention of 5-fold cross validation.

    model_evaluate.m: Plotting the ROC curves based on the 5-fold predicted results.

  2) loocv

     main.m: The main program for the implemention of leave-one-out cross validation..

     Data_fusion_SNF.m: Integrate the similarity matrix of lncRNA/miRNA.

     SNFHGILMA_LOOCV.m: The subprogram for the implemention of leave-one-out cross validation.

     Plot_roc_curve.m: Plotting the ROC curves based on the loocv predicted results.
