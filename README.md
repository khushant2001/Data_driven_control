# 1. Using Singular Value Decomposition to recreate the fluid flow over an airfoil!

In order to understand the significance of Singular Value Decomposition, especially in the field of compression, a big dataset of flow velocities(found ![here](http://deepblue.lib.umich.edu/data/collections/kk91fk98z) is subjected to the algorithm. Specifically, the dataset for airfoil with angle of attack of 25 degrees and frequency of 0.05 hz is used. The dataset itself (for both X and Y velocity) comes in a 3-D matrix, where the first two dimensions represent the velocity for individual snapshots in time. Accordingly, the third dimension contains the number of snapshots, which are 400 in total. Before performing SVD on the dataset, pre-processing is done in the following format:

- The 2-D grid of velocities is shaped into a single column vector.
- Individual snapshots are then horizontally stacked as column vectors. Accordingly, the final product is a big 2-D matrix where columns are snapshots of velocities at consecutive times.
- The mean of every column is subtracted, serving the purpose of normalizing the data.

Finally, the data is ready to be subjected to an economy SVD and the following matrices are derived. 
$$
U*Z*V' = svd(data)
$$
