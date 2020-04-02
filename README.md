# DiscriminantMotifs
Reproducibility of the results from ''Discriminatory Motifs of Complex Networks''  and a few more things that have been analysed but that are not in the paper

Architecture :
Folders :
    -> Processing/Networks :
        Networks on the following format
            # one or several lines
            # that give indication about
            # the network
            !n:number_of_nodes
            !m:number_of_egdes
            v_src v_tgt
            ...
      v_src v_tgt is 1 edge whit v_src source node, and v_tgt target node.
      Node's labels are integers, from 0 to n. The letters before the first
      '-' indicate the field of the network (fw:  foodweb, elec: electronic
      circuits, stac: discourse structure, soc: social networks).
    -> Processing/CountMotifs
       decomposition of networks into 3-node and 4-node motifs. Files are
       labelled such that : {nwk-id}_motif{#k}-counts.tab, where
       nwk-id is the name of the file containing the decomposed network in
       Processing/Network, #k =3 or 4 (3-node opr 4-node decomposition).
       The files have the same format than when using motif from SNAP, but
       they have been generated with Acc-Motif, and then transformed onto
       this format.

MatNetworks: The structure that characterise the networks, generated by the
read_graph script (see the script for further info)

Current : (Only the headlines, more formal information can be found as
comments in the listed files)

    -> read_graph : transforms the networks from the characteristics in
       Processing folder to the mat structure in MatNetwork, which will be
       the format used in all the remaining matlab files.

    -> stats : to generate Table 1.

    -> run_pca_score : code to compute scores in Tables 2 and 3 and display
       all figures from Section3.1

    -> discover_discriminant_motifs : code to generate Fig3

    -> choose_motifs : not used in the paper. Code used to find the number
       of discriminant motifs to keep.

    -> run_motif_score : code to generate scores in Table4.

    -> motif_dispersion : code to generate Fig4.

    -> compare_to_known_motifs : to observe the distribution of networks
       according to a set of classical motifs such as bifan, feedforward
       loop, etc. This code has been used to provide Figures from Sec4.

    -> print_graph : in case one wants to check which graph is associated
       with a given motif identifier (see Sec3 of Supplementary).

Other scripts and functions are used in these main scripts.
