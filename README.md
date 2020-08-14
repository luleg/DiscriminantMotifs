# DiscriminantMotifs
Reproducibility of the results from ''A Simple Embedding for Classifying Networks with a few Graphlets''

## Generation of Tables and Figures

### Generate Table I
  run stats.m

### Generate Table II
#### Our method :
  run run_pca.m
#### Gl2Vec :
  run run_gl2vec.m
#### Graph2Vec :
  run run_graph2vec.m

### Generate Figure I
#### Gamma-score :
  run run_GammaAnalysis.m
#### RF-score :
  run run_RFSelection.m

### Generate Table II
#### Our method :
  run run_feature_select.m
#### RF feature selection :
  run run_RFSelection.m

### Generate Figure II and Confusion Matrices from Table IV
#### Figures Affinity Matrix and Threshold Sparsification / Confusion Matrix to the Left :
  run afty_threshold.m
#### Figure Closest Neighbour Sparsification / Confusion Matrix to the Right :
  run afty_knn.m

## Architecture :

### Folder Stock :
  Matlab files, scripts and functions used by the main run_* scripts
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
#### MatNetworks :
  Each mat file contains a struct Pbm that contains information about the network:
  
    -> Pbm.entete : textual information (website, preprocessing, etc.)
    -> Pbm.nb_nodes/nb_edges : number of nodes/edges (a bidirected edges counts for two edges)
    -> Pbm.motif3 : a matrix 13x2. Pbm.motif3(k,1) : id of 3-node kth motif
                                 Pbm.motif3(k,2) : occurrence number of motif k in the network
    -> Pbm.motif4 : same for 4-node motif
    -> Pbm.edges : a matrix Pbm.nb_edges x 2 where (Pbm.edges(i,1), Pbm.edges(i,2)) = (v_srci,v_tgti)

### Gl2Vec :
#### MotifCounts
  Outputs of the java code from https://github.com/kuntu/JGraphlet-JMotif for our benchmarks.
#### EmbMat
  The SRPs of each networks in Matlab files (generated using convert2Mat.m)
### Graph2Vec :
  A Python code using NetworkX and karate-club Benchmarks to generate the embeddings using graph2vec, for deep of WL-kernel from 1 to 4 and embedding size from 2 to 512.
#### Embedding
  Outputs of the Python code
#### EmbMat
  The corresponding Matlab files (generated using convert_csv2mat.m)

### FeatureSelectionRF
  A Python Code to obtain the average Gini importance of each graphlets by training a forest of 100 trees (using scikit-learn).

### Clusterix :
  Matlab files to run the unsupervised clustering algorithm used in Section VI.
  (See https://pdfs.semanticscholar.org/6235/cf4b551f768fa793ed759de75f2a01475e77.pdf for the algorithm description).
