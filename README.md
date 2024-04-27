Under this project, data driven numerical methods are used to perform system identification and model reduction. The data sets used are subjected to basic linear algebra techniques such as Principle Component Analysis (PCA) and Singular Value Decomposition (SVD) for model reduction. This is followed by system id using methods such as Dynamic Mode Decomposition (DMD), Eigenanalysis Realization Algorithm (ERA), and Sparse Identification of Non Linear Dynamics (Sindy) to generate linear models for non-linear systems. Since the main goal is to design controllers for these systems, linear systems allow to construct controllers such as Linear Quadratic Controller (LQR) and Model Predictive Control (MPC).

# 1. Using Singular Value Decomposition to recreate the fluid flow over an airfoil!

In order to understand the significance of Singular Value Decomposition, especially in the field of compression, a big dataset of flow velocities(found [here](http://deepblue.lib.umich.edu/data/collections/kk91fk98z) is subjected to the algorithm. Specifically, the dataset for airfoil with angle of attack of 25 degrees and frequency of 0.05 hz is used. The dataset itself (for both X and Y velocity) comes in a 3-D matrix, where the first two dimensions represent the velocity for individual snapshots in time. Accordingly, the third dimension contains the number of snapshots, which are 400 in total. Before performing SVD on the dataset, pre-processing is done in the following format:

- The 2-D grid of velocities is shaped into a single column vector.
- Individual snapshots are then horizontally stacked as column vectors. Accordingly, the final product is a big 2-D matrix where columns are snapshots of velocities at consecutive times.
- The 2-D grid of X and Y velocity is stacked vertically to make the analysis less repetetive. 
- The mean of every column is subtracted, serving the purpose of normalizing the data.

Finally, the data is ready to be subjected to an economy SVD and the following matrices are derived: U, S, V. The matrix U contains information about the spatial modes while the product of S and V' contains information about the temporal ampltiude of the corresponding modes in U. As expected, SVD figures out which pattern in the data summarizes it the best. A graphic of the eigenvalues (square of the singular values of the data) is show below. 

![image](https://github.com/khushant2001/Data_driven_control/assets/70731991/64407632-568d-4aed-93bb-6b4e699da9dd)


The first 6 spatial modes which are basically the first 6 column vectors of U are shown below as well. Since the big data matrix contains both x and y velocities, it is important to remember that the top half of the matrix contains the x velocity and bottom half contains the y velocity.   
![image](https://github.com/khushant2001/Data_driven_control/assets/70731991/bab8c3ad-800c-4c7e-8a6e-0229efa5a7b3)

![image](https://github.com/khushant2001/Data_driven_control/assets/70731991/4050b88c-b6e6-4b1b-8012-da7c305489ad)

The corresponding temporal amplitude for these modes is also shown below. As seen, the first 2 dominant modes follow a cyclic wave which are at an offset of a quarter wave. 

![image](https://github.com/khushant2001/Data_driven_control/assets/70731991/9bf44024-9b47-4488-9e95-893772ad9b4f)

Finally, the construction of the flows is performed by choosing a certain number of columns (<400). The mean subtracted from individual columns is added back to the columns selected for reconstruction and the columns are also shaped back to the grid format. By choosing which snapshot to view, a contour plot can be generated, as shown below for first snapshot or t = 1 sec. 

![image](https://github.com/khushant2001/Data_driven_control/assets/70731991/ab65bed5-df53-49bd-9f9c-f675f4b9e0e7)

![image](https://github.com/khushant2001/Data_driven_control/assets/70731991/96593a46-84d4-481b-a1ed-6aa8f1f06463)
A [video](FluidFlowReconstruction_SVD.mp4) of the reconstruction of the flows can be found in the files attached in the repository. 

# 2. Using Dynamic Mode Decomposition to recreate the fluid flow over an airfoil!

Dynamic mode decomposition tries to find the best fit "linear operator" which can predict the next state of the system. Consider 'X' to be the state of the system at time t=k and "X'" to be the state at k+1, then DMD finds the A that can satsify the equations listed below. By taking the inverse of X', A can be determined. 

$$
X' = A * X
$$

$$
A = X' * X^{-1}
$$

Since the matrix X is not always square, Singular Value Decomposition from above is used to find the pseudo-inverse. Equations below show expansion of the matrix X based on SVD and its inverse in terms of the individual matrix. Another important characteristic about DMD is that it finds a matrix called A tilda which is much smaller in size as compared to the original A but has the same eigenvalues. This decreases the computation cost by a lot.  

$$
X = U * \Sigma * V*
$$

$$
X' = A * U * \Sigma * V*
$$

$$
A\sim = U * X' * V * \Sigma^{-1} = U* * A * U
$$

Once A tilda is know, its eigenvalues (W) can be calculated which are same as A and then the eigenvectors of A which are labelled as phi can be determined through back calculation using lambda.  

$$
A\sim * \lambda = \lambda * W
$$

$$
\phi = X' * V * \Sigma^{-1} * \lambda
$$

To perform a DMD on the same flow field, the velocities are stacked in a similar format as before to make the X matrix and the first column is deleted to make X'. The airfoil data used for this has an angle of attack of 25 degrees with freq of 0.05 Hz. The original flow field for the flow at t = 1 looks like below: 
![image](https://github.com/khushant2001/Data_driven_control/assets/70731991/13852d03-5605-4446-96e9-d64bfc47e1b2)

Subjecting to DMD analysis, the amplitudes of the modes are derived and the plot is listed below: 
![image](https://github.com/khushant2001/Data_driven_control/assets/70731991/609e32d1-e82c-40ab-b6bc-28cc80ca54e5)

Also, the contour plot for the first 4 modes along with their frequency is given below: 
![image](https://github.com/khushant2001/Data_driven_control/assets/70731991/9cdd2b8f-059f-40f5-a512-b22349c429e8)

Finally, the reconstruction of the fluid flow is performed using the DMD Algorithm and the following the comparison: 
![image](https://github.com/khushant2001/Data_driven_control/assets/70731991/ff2f27e4-d759-46b6-b774-48ff3b24ab54)

[Watch the video](FluidFlowReconstruction_DMD.mp4)


# 3. Using Sparse Identification of Non-Linear Dynamics (Sindy) to approximate the differential equations for the flow fields.

![image](https://github.com/khushant2001/Data_driven_control/assets/70731991/30c31eea-a137-415e-acb4-fbd71c2d388c)

![image](https://github.com/khushant2001/Data_driven_control/assets/70731991/46e38a60-2446-461b-824d-778bac87712c)

# 4. Eigensystem Realization Algorithm (ERA)
ERA is a data based analysis of a system which is capable of determining the discrete state space matrices: A, B, & C using the impulse response. Accordingly, the the system can be described as follows: 

$$
x_\text{k+1} = A_d x_k + B_d u_k 
$$

$$
y_k = C_d x_k
$$

$$
u_k = \delta(t)
$$

By writing out the response of the system with initial conditions to be 0, the output can be characteristed by 

$$
y_k = C_d A_d ^ {k-1}B_d
$$

Following this the Henkel matrices -H and H' - are generated and the SVD is taken to convert the model from high to low dimensional model. A, B, and C matrices are generated using the following:

$$
SVD(H) = U \sigma V^{*}
$$

$$
A = \sigma^{-.5} U^{T} H' V \sigma^{.5}
$$

$$
B = \sigma^{.5} U^{T} H 
$$

$$
C = H V \sigma^{.5}
$$

To understand this approximation, ERA is tried on the transfer function decribing the Atomic Force Microscope. The corresponding transfer function is given below: 

$$
G = \frac{kw_2^2w_3^2w_5^2(s^2+2\zeta_1w_1s+w1^2)(s^2+2\zeta_4w_4s+w_4^2)exp(-s\tau)}{w_1^2w_4^2(s^2 + 2\zeta_2w_2s + w_2^2)(s^2 + 2\zeta_3w_3s + w_3^2)(s^2 + 2\zeta_5w_5s + w_5^2)}
$$

The impulse response of the system is determined (also shown below) and hankel matrices are generated. 

![image](https://github.com/khushant2001/Data_driven_control/assets/70731991/c1d06d41-6d2b-480b-a903-6dc53d15a4c9)

After determining the discrete A,B,C matrices, the bode plots are generated from the ERA esimtation and actual transfer function. The results are shown below. 
![image](https://github.com/khushant2001/Data_driven_control/assets/70731991/22efa4f7-96e3-4e44-9ef2-29b5e0b53331)

To see the effects of the addition of the white to the system, random normal gaussian noise is added to the impulse response and subjected to the ERA estimation. Surprisingly, the noise does not decrease the accuracy of the estimation; only does at high frequency. This makes sense because by taking SVD of the hankel matrix, the low energy pod modes are neglected which contain the noise. The results of the addition of white noise to the system are shown below:

![image](https://github.com/khushant2001/Data_driven_control/assets/70731991/68001cd2-2a6b-4af2-b1d0-5cc125229b1e)

![image](https://github.com/khushant2001/Data_driven_control/assets/70731991/3c3d6f50-45ea-400f-9e70-9dfb28aecc0c)

# 5. Using Sindy to design a MPC for a Flight Controller

The dataset used can be found here: [here](https://c3.ndc.nasa.gov/dashlink/resources/294/)
