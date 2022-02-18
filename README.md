**Counting Objects by Diffused Index:geometry-free and training-free approach**

Copyright 2021, All Rights Reserved Code by Mengyi Tang, Maryam Yashtini, Sung Ha Kang for Paper, "Counting Objects by Diffused Index: geometry-free and training-free approach" (https://arxiv.org/abs/2110.08365)

main_CODI_M.m    

     -Counting Objects by Diffused Index: Scalar Seeds


main_CODI_S.m    

     -Counting Objects by Diffused Index: Vector Seeds


normalize(x)

    -A helper method to return the normalized data of input in [0,255]
    -inputs:
     x            : Input matrix of size m, n
    -output :
     u            : Output matrix of size m, n, values in [0,255]


add_frame(phi,k, value)

    -A helper method to add a frame of some value and width k into a matrix
    -inputs:
     phi          : a gray-scale image of size M * N
     k            : the width of the frame
     value        : the value of the boundary
    -output:
     framed_phi   : a grapy-scale image of size (M+k) * (N+k)


Treshold2BW_up (x,t)

    -A helper method to return a binary matrix of value 0 and 255 of the
     input matrix
    -inputs:
     x            : a matrix of size M * N
     t            : threshold
    -output:
     u            : a matrix of size M * N, with value 0 where x(i,j) < t
                                            with value 255 o.w.



AlgIEdge(SeedIm,EdgeIm,X_D,mu,theta, eta, perc)

    -A helper function about the Edge-weighted harmonic variational optimization model
     proposed in the paper to produce a diffused image in a single
     dimension
    -inputs:
     SeedIm       : Seed Image
     EdgeIm       : Mask Image of the image that $g$ applied onto
     mu           : multiplier on the constraint V−U=0.
     lambda       : λ, the Lagrange multiplier associated with the linear constraint V−U=0.
     eta          : η, fidelity parameter
    -output :
     DiffiusedIm  : Diffused image
     objval       : Energy Function
     Itercount    : Iteration Number
     CPUcount     : CPU time



PlotClusterinResult(D, IDX, dir1, dir2)

    -A helper method to plot the clustering results from 4 dimensional data
     into a subspace contructed by dir1 and dir2 dimension
    -inputs:
     D            : a matrix of size m*n, m data points in n dimensions
     IDX          : a list of cluster index of each point in D
     dir1, dir2   : the index of the subspace, default values are 1 and 2
    -output:
     Figure       : projection of clustering results from high dimension to
                    low dimension
