# Bjorck-Duff

This repository contains the code scripts and other related matrials for implementing Bjorck-Duff's algorithm [1]. This algorithm solves the linear least-squares problem $\min_x || b - Ax ||_2$ and can be viewed as an adaptation of Peters and Wilkinson's algorithm [2] to sparse matrices. We discuss how we ensure numerical stability while preserving much of the sparsity of the original system.

Furthermore, we have implemented the algorithm in MATLAB using test sparse matrices from the [SuiteSparse matrix collection](http://sparse.tamu.edu/).

---------

[1] Björck, Åke, and Iain S. Duff. "A direct method for the solution of sparse linear least squares problems." _Linear Algebra and its Applications_ 34 (1980): 43-67.

[2] Peters, Gwen, and James Hardy Wilkinson. "The least squares problem and pseudo-inverses." _The Computer Journal_ 13, no. 3 (1970): 309-316.

---------

(c) Praveen Kumar, Rohan Karthikeyan, Roudranil Das, Saikat Bera.

__*This code is being released with an objective to enhance general understanding of the research paper and its underlying concepts. Unattributed publishing and uploading of these contents with or without modification in the Internet shall not be encouraged.*__
