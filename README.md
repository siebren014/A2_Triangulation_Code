# A2_Triangulation_Code
GEO1016_Assignment_2_Group_2 description

README triangulation_method.cpp structure 

Line 33
//----requirements for input data
// we need at least 8 pairs for the 8-point algorithm
Bool isvalid(const std::vector<Vector2D> &points_0, const std::vector<Vector2D> &points_1)

Line 45
//---- translation and scaling matrix
std::pair<Matrix33, bool>Transform_mat_normalized(
        const std::vector<Vector2D>& points)
Line 93
//----Normalized eight-point algorithm
std::vector<std::vector<Vector2D>> normalize_points( const std::vector<Vector2D> &points_0,const std::vector<Vector2D> &points_1) {

Line 160
//----- determination of the 3D points constructed from both images and their possible R and t
std::vector<Vector3D> Points(const Matrix33 &K, const std::vector<Vector2D> &points_0,const std::vector<Vector2D> &points_1, const Matrix &R, const Vector &t)

Line 210
//----- counter for the different point_options for the different R and t
std::vector<std::vector<Vector3D>> point_options (const Matrix33 &K, const std::vector<Vector2D> &points_0,const std::vector<Vector2D> &points_1, const Matrix &R1, const Vector &t1, const Matrix &R2, const Vector &t2)

Line 226
//----- Reading in the different points and construct the F and E matrix
bool Triangulation::triangulation

Line 333
//----- searching for the best combination of R and t
auto point_candidates = point_options(K, points_0, points_1, R1, t1, R2, t2);

Line 336
// option_score keeps track of the combination of the cases in which ppoint.z>0
std::vector<int> options_score;

Line 375
// maxelementindex checks the options_score in the beginning and the end
int maxElementIndex = (std::max_element(options_score.begin(),options_score.end()) - options_score.begin());
points_3d = point_candidates[maxElementIndex];

Line 401
//---- evaluation of the results    
