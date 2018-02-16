# T-LoCoH_Algorithm

The algorithm, as presented in Movement Ecology (https://doi.org/10.1186/s40462-017-0110-4), carries out a grid-based search of k and s values in an attempt to optimize the T-LoCoH parameter selection process. This particular approach values consistency in home range delineations (as measured through cross-validation of alternative training-testing sets) over alternatives, such as avoiding a swiss-cheese appearance.

The two files here consist of the original algorithm as included in the Supplementary Materials of the paper above and another version that divides the procedure into two component parts (for the sake of replicability): the training/testing split process and the actual grid-based search. See the code itself for more details about the arguments accepted by each.
