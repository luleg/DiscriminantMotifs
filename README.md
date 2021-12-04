# Discriminate Networks Using a Few Graphlets
Reproducibility of the results from ''A Simple Embedding for Classifying Networks with a few Graphlets''

## Generation of Tables and Figures

### Generate Table I
  run ``stats.m``

### Generate Table II
#### PCA :
  run ``run_pca.m``
  
  To change the normalisation, modify global variable ``type_norm`` (line 43)
#### LDA :
  run ``run_lda.m``
  
  To change the normalisation, modify global variable ``type_norm`` (line 43)

### Generate Table III
  run ``run_pca.m``
  
  To change the kind of k-node graphlets to be used, modify ``motif3/4/5`` (lines 16-18)

### Generate Table IV
  run ``run_graph2vec.m``
  
  To change the size of the embeddings, modify the global variable ``t_emb`` (line 21)
  
    Precomputed embedding sizes are 2, 4, 8, 16, 32, 64, 128, 256, 512.
  
  To change the depth of the WL kernel, modify the global variable ``t_wl`` (line 21)
  
    Precomputed WL depths are 1, 2, 3

### Generate Table V
  Go in the GCNs folder and create the dataset : ``python3 create_graphtest.py``
#### GCNs
  Train and run the neural network :``python3 run_GCN.py -pooling max/mean/sum -hidd 4/8/12 -nbLay 1/2/3``
  
      You may also change the variable ``nbItes`` (line 15), so that ``nbItes = 10``
  
  Obtain scores : ``python3 analysis.py GCN_poolmax/mean/sum_nbLay1/2/3_hidd4/8/12``
#### RGCNs
  Train and run the neural network : ``python3 run_RGCN.py -pooling max/mean/sum -hidd 4/8/12 -nbLay 1/2/3``
  
      You may also change the variable ``nbItes`` (line 15), so that ``nbItes = 10``
  
  Obtain scores : ``python3 analysis.py GCN_poolmax/mean/sum_nbLay1/2/3_hidd4/8/12``


### Generate Table VI
#### Our method :
  run ``run_pca.m``
#### Gl2Vec :
  run ``run_gl2vec.m``
#### Graph2Vec :
  run ``run_graph2vec.m``
#### GCNs (takes a long time)
  Train and run the neural network : ``python3 run_GCN.py``
  
      Change the variable ``nbItes = 50`` (line 15)
      
  Obtain scores : ``python3 analysis.py GCN_poolmean_nbLay2_hidd12``
#### RGCNs (takes a long time)
  Train and run the neural network : ``python3 run_RGCN.py``
  
      Change the variable ``nbItes = 50`` (line 15)
  
  Obtain scores : ``python3 analysis.py RGCN_poolmean_nbLay2_hidd8``


### Generate Table VII
  run ``run_pca_wiki.m``

### Generate Figures VI, VII, VIII
  run ``run_Wiki_figures.m``

### Generate Figure IX
#### Gamma-score :
  run ``run_GammaAnalysis.m``
#### RF-score :
  run ``run_RFAnalysis.m``

### Generate Table VIII
#### Our method :
  run ``run_feature_select.m``
#### RF feature selection :
  run ``run_RFAnalysis.m``

### Generate Figure X and Confusion Matrices from Table X
#### Figures Affinity Matrix and Threshold Sparsification / Confusion Matrix to the Left :
  run ``afty_threshold.m``
#### Figure Closest Neighbour Sparsification / Confusion Matrix to the Right :
  run ``afty_knn.m``

## Architecture :

### Folder Stock :
  Matlab files, scripts and functions used by the main ``run_*`` scripts
#### Processing/Networks :
  Networks used in the tests, on the following format :

    # one or several lines
    # that give indication about
    # the network
    !n:number_of_nodes
    !m:number_of_egdes
    v_src1 v_tgt1
    v_src2 v_tgt2
    ...
#### Processing/Wiki_Graph_remove_40_0_500 :
  Wikipedia Networks, on the following format :

    # Article : Name_of_the_Wikipedia_Article
    # FROM date_of_the_first_version_of_interest TO date_of_the_last_version_of_interest
    number_of_nodes number_of_edges
    v_src1 v_tgt1
    v_src2 v_tgt2
    ...

#### MatNetworks :
  Each mat file contains a struct Pbm that contains information about the network:
  
    -> Pbm.entete : textual information (website, preprocessing, etc.)
    -> Pbm.nb_nodes/nb_edges : number of nodes/edges (a bidirected edges counts for two edges)
    -> Pbm.motif3 : a matrix 13x2. Pbm.motif3(k,1) : id of 3-node kth motif
                                   Pbm.motif3(k,2) : occurrence number of motif k in the network
    -> Pbm.motif4 : same for 4-node motif
    -> Pbm.motif5 : same for 5-node motif (does not exist for all networks)
    -> Pbm.edges : a matrix Pbm.nb_edges x 2 where (Pbm.edges(i,1), Pbm.edges(i,2)) = (v_srci,v_tgti)

### Gl2Vec :
#### MotifCounts
  Outputs of the java code from https://github.com/kuntu/JGraphlet-JMotif for our benchmarks.
#### EmbMat
  The SRPs of each networks in Matlab files (generated using ``convert2Mat.m``)
### Graph2Vec :
  A Python code using NetworkX and karate-club Benchmarks to generate the embeddings using graph2vec, for deep of WL-kernel from 1 to 4 and embedding size from 2 to 512.
#### Embedding
  Outputs of the Python code
#### EmbMat
  The corresponding Matlab files (generated using ``convert_csv2mat.m``)

### FeatureSelectionRF
  A Python Code to obtain the average Gini importance of each graphlets by training a forest of 100 trees (using scikit-learn).

### Clusterix :
  Matlab files to run the unsupervised clustering algorithm used in Section V.
  (See https://pdfs.semanticscholar.org/6235/cf4b551f768fa793ed759de75f2a01475e77.pdf for the algorithm description).
