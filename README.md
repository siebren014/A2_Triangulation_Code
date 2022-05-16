# A2_Triangulation_Code
GEO1016_Assignment_2_Group_2 description

step 1 estimate fundamental matrix F, from given corresponding image points using the normalized 8point algorithm (20%)
a. normalization                                                (5%)
b. linear solution (based on SVD), find closes rank 2 matrix    (5%)
c. constraint enforcement based on SVD                          (5%)
d. denormalization                                              (5%)

step 2. recover relative pose (i.e. R and t) from the fundamental matrix (20%)
a. find the 4 candidate relative poses (based on SVD)           (10%)
b. determine the correct relative pose                          (10%)
(in front of the camera)

step 3. determine the 3D coordinates for all corresponding image points (20%)
a. compute the projection matrix from K, R and t                (5%)
b. compute the 3D point using the linear method (based on SVD)  (5%)
b+ (optional) non linear least squares refinement of the 3D point computed from the linear method (10% bonus)
c. triangulate all corresponding image points                   (10%)

step 4 evaluation (10%)

step 5 report (30%)
a. description methodology                                      (5%)
b. demonstration of your result                                 (5%)
c. evaluation of your results (accuracy and a way to improve)   (5%)
d. discussion on how to obtain the intrinsic parameters         (5%)
e. reflection of another group                                  (5%)
f. short description of task division                           (5%)
+ a file of the reconstructed 3D points (in .xyz format)
