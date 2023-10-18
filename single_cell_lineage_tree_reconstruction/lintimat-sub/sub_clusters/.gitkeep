# *src/Data/* containts the dataset  

The following are the required datasets  

    1. Homing barcode mutations  
    2. Mitochondrial mutations  
    3. RNA data (.h5ad)
    4. Cell-type annotations of the cell

# Run

Run the `src/lintimat/code/Lintimat.ipynb` notebook.

## *src/lintimat/code/Lintimat.ipynb*  

contains the code for the compleate generation of the output files in the *src/lintimat/output/* folder. While the processed files required by the lintimat get stored in the *src/lintimat/input/*  folder [ref lintimat](https://github.com/jessica1338/LinTIMaT/wiki/Input-File)

## The steps in Lintimat.ipynb

    1. Pre-process homing and Mitochondrial data.
    2. Pre-process RNA data.
    3. Create cell type mapping color file.
    4. Create dataset for lintimat.
    5. runs the lintimat.
    6. Generation and visualization of the non-binary tree.  

Details in the methods session.

# For larger datasets

The following are the required datasets.

    1. Homing barcode mutations  
    2. Mitochondrial mutations  
    3. RNA data (.h5ad)
    4. Cell-type annotations of the cell

## 1. src/lintimat-sub/code/Sub trees.ipynb

The dataset is pre-processed and divided into sub-clusters based on the mutation profiles and lintimat trees are generated for each sub-cluster (sub-trees). The above steps are followed in `Sub trees.ipynb` notebook.  
The processed files are saved in the `src/lintimat-sub/input/`, and the `src/lintimat-sub/sub_clusters/` will store the lintimat results for each sub cluster.

## 2. src/lintimat-sub/code/Compleate tree.ipynb 

Further, a backbone tree is generated based on each sub-cluster's binarized average mutation profile using lintimat. Then, each sub tree is attached to the end of the backbone tree to get the compleate tree.