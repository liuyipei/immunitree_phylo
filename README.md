Immunitree is implemented in matlab
============
**Contact:**
 - liu.yi.pei@gmail.com
 - jonilaserson@gmail.com
**Reference:**
 - http://www.plospathogens.org/article/info%3Adoi%2F10.1371%2Fjournal.ppat.1003754
 - https://stacks.stanford.edu/file/druid:xp796hy4748/thesis-augmented.pdf

After downloading the entire repository from git
==========
 - Examine `Fire/sample_driver_script.m`
 - Update `Fire/sample_driver_script.m` with correct paths to the location of the `immunitree_phylo directory`
 - Update `Fire/sample_driver_script` with correct paths to the locations of the input files
 -  Make sure that the chain type parameter (0 for heavy, 1 for light) is correct provided to the call to `pipeline_for_hiv.m`
 - Examine the `.../Fire/output/* folder for output`
 - `*.txt` files: descriptions of the trees for parsing (There are 3 txt files per tree)
  - The three txt files corresponds to three perspectives of the tree: overall, per-node, and, per-read
  - There are also png files to visualize the trees. The colorings correspond to the number of mutations from the assigned node
