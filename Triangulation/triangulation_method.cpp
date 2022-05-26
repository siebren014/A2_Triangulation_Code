/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>
#include <cmath>


using namespace easy3d;

//----requirements for input data
// we need at least 8 pairs for the 8-point algorithm
bool isvalid(const std::vector<Vector2D> &points_0, const std::vector<Vector2D> &points_1){

    if (points_0.size() >= 8 && points_0.size() == points_1.size()){
        return true;
    }
    else{
        return false;
    }
}

//---- translation and scaling matrix 
std::pair<Matrix33, bool>Transform_mat_normalized(
        const std::vector<Vector2D>& points)
{
    //construction of the T matrix 
    Matrix33 T;
    bool T_valid = false;
    //getting the centre of the points
    double sumX = 0;  
    double sumY = 0;
    const double N = static_cast<double>(points.size());

    for (const auto& p : points)
    {
        sumX += p.x();
        sumY += p.y();
    }
    if (sumX < 1e-8 || sumY < 1e-8)  
    {
        LOG(ERROR) << "error input\n";
        return std::make_pair(T, T_valid);  /* T_valid is false, no further process steps*/
    }

    //centre of the image 
    double tx = sumX / N;
    double ty = sumY / N;

    //scaling factor 
    double sum_dist = 0;
    for (const auto& p : points)
        sum_dist += sqrt((p.x() - tx) * (p.x() - tx) + (p.y() - ty) * (p.y() - ty));
    double avg_dc = sum_dist / N;
    if (avg_dc < 1e-8)
    {
        LOG(ERROR) << "no correct average of the points\n";
        return std::make_pair(T, T_valid);   /* T_valid is false, no further process steps*/
    }
    double s = sqrt(2) / avg_dc;
    //construction of the transformation matrix 
    T.set_row(0, { s, 0, -s * tx });
    T.set_row(1, { 0, s, -s * ty });
    T.set_row(2, { 0, 0, 1 });
    T_valid = true;  /* T is successfully constructed */

    return std::make_pair(T, T_valid);
}


//----Normalized eight-point algorithm
std::vector<std::vector<Vector2D>> normalize_points( const std::vector<Vector2D> &points_0,const std::vector<Vector2D> &points_1) {
    std::vector<std::vector<Vector2D>> normalize_results;
    std::vector<Vector2D> points_0_normalized;
    std::vector<Vector2D> points_1_normalized;
    Vector2D sum_points0;
    Vector2D sum_points1;
    Vector2D sum_points0norm;
    Vector2D sum_points1norm;
    double dis0 = 0;
    double dis1 = 0;
    double dis0_norm = 0;
    double dis1_norm = 0;

    for (const auto &p1: points_0) {
        sum_points0 += p1;
    }
    for (const auto &p1: points_1) {
        sum_points1 += p1;
    }

    //mean position calculation of the both camera standpoints, translation
    Vector2D mean_points0 = sum_points0 / points_0.size();
    Vector2D mean_points1 = sum_points1 / points_1.size();

    std::vector<Vector2D> tr_points0;
    std::vector<Vector2D> tr_points1;
    
    // all the points translated to their image centre 
    for (const auto &p1:points_0){
        tr_points0.emplace_back(p1 - mean_points0);
    }
    for (const auto &p1:points_1){
        tr_points1.emplace_back(p1 - mean_points1);
    }

    //eventually scaling of the points 
    //determination of their distance from the centre of their image
    for (const auto &p1: tr_points0) {
        dis0 += distance(p1, Vector2D(0,0));
    }
    auto avg_dis0 = dis0 / points_0.size();

    for (const auto &p1: tr_points1) {
        dis1 += distance(p1, Vector2D(0,0));
    }
    auto avg_dis1 = dis1 / points_1.size();

    //average distance divided by the sqrt of 2, for normalization
    auto norm_factor = avg_dis0 / sqrt(2)  ;
    for (auto &p1: tr_points0) {
        points_0_normalized.emplace_back(p1 / norm_factor);
        sum_points0norm += p1 / norm_factor;
    }
    auto norm1_factor =  avg_dis1/ sqrt(2);
    for (auto &p1: tr_points1) {
        points_1_normalized.emplace_back(p1 / norm1_factor);
        sum_points1norm += p1 / norm1_factor;
    }

    normalize_results.emplace_back(points_0_normalized);
    normalize_results.emplace_back(points_1_normalized);

    return normalize_results;
}


//----- determination of the 3D points constructed from both images and their possible R and t

std::vector<Vector3D> Points(const Matrix33 &K, const std::vector<Vector2D> &points_0,const std::vector<Vector2D> &points_1, const Matrix &R, const Vector &t){
    Matrix r_t = Matrix(3,4);
    r_t.set_column(0,{R.get_column(0)});
    r_t.set_column(1,{R.get_column(1)});
    r_t.set_column(2,{R.get_column(2)});
    r_t.set_column(3,{t});

//    std::cout<<"r_t"<<r_t<<std::endl;

    Matrix identit = Matrix(3,4);
    identit.set_row(0,{1,0,0,0});
    identit.set_row(1,{0,1,0,0});
    identit.set_row(2,{0,0,1,0});

    // M=K[IO] for camera 0 and M=K[Rt] for camera 1
    auto M_prime = K*r_t;
    auto M = K * identit;

    //triangulation method to determine the specific 3D points
    std::vector<Vector3D> points_3d;
    for (int i = 0; i < points_0.size(); i++){

        auto x = points_0[i][0];
        auto y = points_0[i][1];
        auto x_prime = points_1[i][0];
        auto y_prime = points_1[i][1];

        Matrix44 A;
        A.set_row(0,{x*M.get_row(2) - M.get_row(0)});
        A.set_row(1,{y*M.get_row(2) - M.get_row(1)});
        A.set_row(2,{x_prime*M_prime.get_row(2) - M_prime.get_row(0)});
        A.set_row(3,{y_prime*M_prime.get_row(2) - M_prime.get_row(1)});
        
        //SVD of matrix A for 3D points 
        Matrix U = Matrix(A.rows(),A.cols(),0.0);
        Matrix S = Matrix(A.rows(),A.cols(), 0.0);
        Matrix V = Matrix(A.cols(),A.cols(),0.0);

        svd_decompose(A,U,S,V);

        Vector4D p = V.get_column(V.cols() - 1);  
        points_3d.emplace_back();
        points_3d.back() = p.cartesian();  
    }

    return points_3d;
}

    //----- counter for the different point_options for the different R and t 
std::vector<std::vector<Vector3D>> point_options (const Matrix33 &K, const std::vector<Vector2D> &points_0,const std::vector<Vector2D> &points_1, const Matrix &R1, const Vector &t1, const Matrix &R2, const Vector &t2){

    std::vector<std::vector<Vector3D>> pointoptions;
    auto option_1 = Points(K, points_0, points_1, R1, t1);
    auto option_2 = Points(K, points_0, points_1, R1, t2);
    auto option_3 = Points(K, points_0, points_1, R2, t1);
    auto option_4 = Points(K, points_0, points_1, R2, t2);
    pointoptions.emplace_back(option_1);
    pointoptions.emplace_back(option_2);
    pointoptions.emplace_back(option_3);
    pointoptions.emplace_back(option_4);

    return pointoptions;
}

    //----- Reading in the different points and construct the F and E matrix 

bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const {

    if (!isvalid(points_0, points_1)) {
        throw std::invalid_argument("invalid input");
    }

    auto norm_points = normalize_points(points_0, points_1);
    auto points_0_norm = norm_points[0];
    auto points_1_norm = norm_points[1];


    // define W_matrix based on amount of inputpoints

    Matrix W_matrix_normalized(points_0.size(), 9, 0.0);

    // fill W_matrix by traversing through all the points.
    for (int i = 0; i < points_0.size(); i++) {
        Vector2D p1 = points_0_norm[i];
        Vector2D p2 = points_1_norm[i];
        auto u1 = p1[0];
        auto v1 = p1[1];
        auto u1prime = p2[0];
        auto v1prime = p2[1];
        W_matrix_normalized.set_row(i,
                                    {u1 * u1prime, v1 * u1prime, u1prime, u1 * v1prime, v1 * v1prime, v1prime, u1, v1,
                                     1});
    }

    int m = points_0.size();
    int n = 9;

    Matrix U = Matrix(m, m, 0.0);
    Matrix S = Matrix(m, n, 0.0);
    Matrix V = Matrix(n, n, 0.0);

    svd_decompose(Matrix(W_matrix_normalized), U, S, V);
    Vector F = V.get_column(V.cols() - 1);

    Matrix F_mat = Matrix(3, 3, 0.0);
    F_mat[0][0] = F[0];
    F_mat[0][1] = F[1];
    F_mat[0][2] = F[2];
    F_mat[1][0] = F[3];
    F_mat[1][1] = F[4];
    F_mat[1][2] = F[5];
    F_mat[2][0] = F[6];
    F_mat[2][1] = F[7];
    F_mat[2][2] = F[8];

    Matrix U_mat = Matrix(F_mat.rows(), F_mat.rows(), 0.0);
    Matrix S_mat = Matrix(F_mat.rows(), F_mat.cols(), 0.0);
    Matrix V_mat = Matrix(F_mat.cols(), F_mat.cols(), 0.0);

    svd_decompose(F_mat, U_mat, S_mat, V_mat);

    //set the rank to 2 with constraint
    S_mat[2][2] = 0;

    //compute the new apprioximated F
    Matrix33 F_bestrank = (U_mat * S_mat * V_mat.transpose());

    auto transform_0 = Transform_mat_normalized(points_0);
    auto transform_1 = Transform_mat_normalized(points_1);

    Matrix33 T = transform_0.first;
    Matrix33 T_prime = transform_1.first;

    Matrix33 F_denorm = T_prime.transpose() * F_bestrank * T;
    auto F_scaled = F_denorm * 1.0 / F_denorm(2, 2);

    // setting the values of the intrinsic parameters
    Matrix33 K;
    K.set_row(0, {fx, 0, cx});
    K.set_row(1, {0, fy, cy});
    K.set_row(2, {0, 0, 1});

    //from fundamental F matrix to essential E Matrix
    Matrix E = transpose(K) * F_scaled * K;

    Matrix U_E = Matrix(E.rows(), E.rows(), 0.0);
    Matrix V_E = Matrix(E.rows(), E.cols(), 0.0);
    Matrix S_E = Matrix(E.cols(), E.cols(), 0.0);

    svd_decompose(Matrix(E), U_E, S_E, V_E);

    //determination of the essential matrix the 4 potential relative pose elements ( 2potential R and 2 potential t).
    
    Matrix W_E = Matrix(3, 3);
    W_E.set_row(0, {0, -1, 0});
    W_E.set_row(1, {1, 0, 0});
    W_E.set_row(2, {0, 0, 1});

    auto R1 = determinant(U_E * W_E * transpose(V_E)) * U_E * W_E * transpose(V_E);
    auto R2 = determinant(U_E * transpose(W_E) * transpose(V_E)) * U_E * transpose(W_E) * transpose(V_E);
    auto t1 = U_E.get_column(U_E.cols() - 1);
    auto t2 = -1 * U_E.get_column(U_E.cols() - 1);

    //----- searching for the best combination of R and t
    auto point_candidates = point_options(K, points_0, points_1, R1, t1, R2, t2);

    // option_score keeps track of the combination of the cases in which ppoint.z>0
    std::vector<int> options_score;
    for (int i = 0; i < point_candidates.size(); i++) {
        int option = 0;
        for (const auto &point: point_candidates[i]) {
            if (point.z()>0 ){
                Matrix33 R_prime;
                Vector3D t_prime;
                switch(i){
                    case 0:
                        R_prime = R1;
                        t_prime = t1;
                        break;
                    case 1:
                        R_prime = R1;
                        t_prime = t2;
                        break;
                    case 2:
                        R_prime = R2;
                        t_prime = t1;
                        break;
                    case 3:
                        R_prime = R2;
                        t_prime = t2;
                        break;

                    default:
                        throw std::invalid_argument( "not 1 unique solution" );
                        ;
                }
                Vector3D ppoint = R_prime * point + t_prime;
                if (ppoint.z()>0 ) {
                    option +=1;
                }
            }
        }
        options_score.emplace_back(option);
    }

    // maxelementindex checks the options_score in the beginning and the end
    int maxElementIndex = (std::max_element(options_score.begin(),options_score.end()) - options_score.begin());
    points_3d = point_candidates[maxElementIndex];
    switch(maxElementIndex) {
        case 0:
            R = R1;
            t = t1;
            break;
        case 1:
            R = R1;
            t = t2;
            break;
        case 2:
            R = R2;
            t = t1;
            break;
        case 3:
            R = R2;
            t = t2;
            break;

        default:
            throw std::invalid_argument( "not 1 unique solution" );
            ;
    }
    
    //---- evaluation of the results 
    //checking the difference between the u,v values in the beginning with the reprojected u,v from the 3D points

    // from the points_options it came forward that the combination of R2 and t1 gave the correct result 
    Matrix rt = Matrix(3,4);
    rt.set_column(0,{R2.get_column(0)});
    rt.set_column(1,{R2.get_column(1)});
    rt.set_column(2,{R2.get_column(2)});
    rt.set_column(3,{t1});

    Matrix identit = Matrix(3,4);
    identit.set_row(0,{1,0,0,0});
    identit.set_row(1,{0,1,0,0});
    identit.set_row(2,{0,0,1,0});

    auto M_prime = K*rt;
    auto M = K * identit;

    //std::cout<<M<<std::endl;

    //calculating the difference between the u,v given and u,v calculated 
    Vector2D diff0;
    Vector2D diff1;

    for (int i = 0; i < points_3d.size(); i++) {
        Vector3D p0_proj = M * points_3d[i].homogeneous();
        Vector2D p0_2d = p0_proj.cartesian();
        Vector3D p1_proj = M_prime * points_3d[i].homogeneous();
        Vector2D p1_2d = p1_proj.cartesian();
        diff0 += p0_2d - points_0[i];
        diff1 += p1_2d - points_1[i];
    }

    std::cout<<"diff0 = "<<diff0/points_3d.size()<<std::endl;
    std::cout<<"diff1 = "<<diff1/points_3d.size()<<std::endl;

    return !points_3d.empty();
}
