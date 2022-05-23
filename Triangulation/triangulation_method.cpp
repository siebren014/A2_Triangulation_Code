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


// we need at least 8 pairs for the 8-point algorithm
bool isvalid(const std::vector<Vector2D> &points_0, const std::vector<Vector2D> &points_1){

    if (points_0.size() >= 8 && points_0.size() == points_1.size()){
        return true;
    }
    else{
        return false;
    }
}

Matrix trans_scale (const std::vector<Vector2D> &points){
    Vector2D sum_points;

    for (const auto &p1: points) {
        sum_points += p1;
    }

    //mean position calculation of the both camera standpoints, translation
    Vector2D mean_points = sum_points / points.size();

    // translation matrix
    auto T = Matrix::identity(3, 3, 1.0);
    T.set(0, 2, -mean_points.x());
    T.set(1, 2, -mean_points.y());

    // scaling part - average distance of translated points from their meanpointsx and meanpointsy
    double dist_sum = 0;
    for (auto &p : points){
        auto homogeneous = p.homogeneous();
        Vector3D translation = T * homogeneous;
        auto dist = distance(translation.cartesian(), mean_points);
        dist_sum += dist;
    }
    double avg_dist = dist_sum / (double)points.size();
    double s = sqrt(2.0) / avg_dist;

    // scaling matrix
    Matrix S_M;
    S_M.set_row(0, (s,0,0));
    S_M.set_row(1, (0,s,0));
    S_M.set_row(2, (0,0,1));

    // transformation matrix from T and S combined. First the translation afterwards scaling
    auto trans_and_scale = S_M * T;
    return trans_and_scale;
};

//Normalized eight-point algorithm
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
    auto tspoints = trans_scale(points_0);
    auto tspoints1 = trans_scale(points_1);

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

    for (const auto &p1:points_0){
        tr_points0.emplace_back(p1 - mean_points0);
    }

    for (const auto &p1:points_1){
        tr_points1.emplace_back(p1 - mean_points1);
    }

    //average distance of the transformed image points from the origin is equal to sqrt2 pixels, scaling
    for (const auto &p1: tr_points0) {
        dis0 += distance(p1, Vector2D(0,0));
    }
    auto avg_dis0 = dis0 / points_0.size();

    for (const auto &p1: tr_points1) {
        dis1 += distance(p1, Vector2D(0,0));
    }
    auto avg_dis1 = dis1 / points_1.size();



    //transfomation matrix T, that translate by the centroid and scale by the scaling factor
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

//    //mean normalized position calculation of the both camera standpoints
//    Vector2D mean_norm_points0 = sum_points0norm / points_0.size();
//    Vector2D mean_norm_points1 = sum_points1norm / points_1.size();
//
//    //average distance of the normalized points ????
//    for (const auto &p1: points_0_normalized) {
//        dis0_norm += distance(p1, mean_norm_points0);
//    }
//    for (const auto &p1: points_1_normalized) {
//        dis1_norm += distance(p1, mean_norm_points1);
//    }
//
//
//    std::cout<<"mean_norm_points0"<<mean_norm_points0<<std::endl;
//    std::cout<<"mean_norm_points1"<<mean_norm_points1<<std::endl;
//    std::cout<<"dis0_norm "<<dis0_norm/160 <<std::endl;
//    std::cout<<"dis1_norm "<<dis1_norm / 160<<std::endl;

    normalize_results.emplace_back(points_0_normalized);
    normalize_results.emplace_back(points_1_normalized);

    return normalize_results;
}



// determine the 4 possible camera positions

std::vector<Vector3D> Points(const Matrix33 &K, const std::vector<Vector2D> &points_0,const std::vector<Vector2D> &points_1, const Matrix &R, const Vector &t){

    Matrix r_t = Matrix(3,4);
    r_t.set_column(0,{R.get_column(0)});
    r_t.set_column(1,{R.get_column(1)});
    r_t.set_column(2,{R.get_column(2)});
    r_t.set_column(3,{t});

    Matrix identit = Matrix(3,4);
    identit.set_row(0,{1,0,0,0});
    identit.set_row(1,{0,1,0,0});
    identit.set_row(2,{0,0,1,0});

    auto M_prime = K*r_t;
    auto M = K * identit;


    std::vector<Vector3D> points_3d;
    // calculating 4d P point
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

        Matrix U = Matrix(A.rows(),A.cols(),0.0);
        Matrix S = Matrix(A.rows(),A.cols(), 0.0);
        Matrix V = Matrix(A.cols(),A.cols(),0.0);

        svd_decompose(A,U,S,V);
        //LEON Vector4D P4 = V.get_column(V.cols() - 1);
        //LEON return P4.cartesian();
        Vector3D p = V.get_column(V.cols() - 1);
        points_3d.push_back(p);
        // NEW Vector4D p = V.get_column(V.cols() - 1);
        // NEW points_4d.push_back(p);
        //points_4d.cartesian();
    }

    //return points_4d.cartesian();
    //return points_3d.pushback(p);
    return points_3d;
}


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


bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const{

    if (!isvalid(points_0, points_1)){
        throw std::invalid_argument( "invalid input" );
    }

    auto norm_points = normalize_points(points_0, points_1);
    auto points_0_norm = norm_points[0];
    auto points_1_norm = norm_points[1];


    //// define W_matrix based on amount of inputpoints

    Matrix W_matrix_normalized(points_0.size(), 9, 0.0);

    //fill W_matrix by traversing through all the points.
    for (int i = 0; i < points_0.size(); i++){
        Vector2D p1 = points_0_norm[i];
        Vector2D p2 = points_1_norm[i];
        auto u1 = p1[0]; auto v1 = p1[1];
        auto u1prime = p2[0]; auto v1prime = p2[1];

//        W_matrix.set_row(i, {points_0[i][0]*points_1[i][0], points_0[i][1]*points_1[i][0], points_1[i][0], points_0[i][0]*points_1[i][1], points_0[i][1]*points_1[i][1], points_1[i][1], points_0[i][0], points_0[i][1], 1});
        W_matrix_normalized.set_row(i, {u1*u1prime, v1*u1prime, u1prime, u1*v1prime, v1*v1prime, v1prime, u1, v1, 1});
    }

    int m = points_0.size();
    int n = 9;

//    Vector F = Vector(n, 0.0);
    Matrix U = Matrix(m,m,0.0);
    Matrix S = Matrix(m,n, 0.0);
    Matrix V = Matrix(n,n,0.0);

    svd_decompose(Matrix (W_matrix_normalized), U, S, V);
    Vector F = V.get_column(V.cols() - 1);

    Matrix F_mat= Matrix(3, 3, 0.0);
    F_mat[0][0]=F[0];
    F_mat[0][1]=F[1];
    F_mat[0][2]=F[2];
    F_mat[1][0]=F[3];
    F_mat[1][1]=F[4];
    F_mat[1][2]=F[5];
    F_mat[2][0]=F[6];
    F_mat[2][1]=F[7];
    F_mat[2][2]=F[8];

    Matrix U_mat = Matrix(m,m,0.0);
    Matrix S_mat = Matrix(m,n, 0.0);
    Matrix V_mat = Matrix(n,n,0.0);

    svd_decompose(F_mat,U_mat, S_mat,V_mat);

    //set the rank to 2 with constraint
    S_mat[2][2]=0;

    //compute the new apprioximated F
    Matrix33 F_bestrank = (U_mat * S_mat * V.transpose());

    auto F_scaled = F_bestrank/ F_bestrank[2][2];

    // setting the values of the intrinsic parameters
    Matrix33 K;
    K.set_row(0,{fx, 0, cx});
    K.set_row(1,{0,fy,cy});
    K.set_row(2,{0,0,1});

    //from fundamental F matrix to essential E Matrix
    Matrix E = transpose(K)*F_scaled*K;

    Matrix U_E = Matrix(E.rows(),E.rows(),0.0);
    Matrix V_E = Matrix(E.rows(),E.cols(),0.0);
    Matrix S_E = Matrix(E.cols(),E.cols(),0.0);

    svd_decompose(Matrix (E), U_E, S_E, V_E);

    //determination of the essential matrix the 4 potential relative pose elements ( 2potential R and 2 potential t).
    Matrix W_E = Matrix(3,3);
    W_E.set_row(0,{0,-1,0});
    W_E.set_row(1,{1,0,0});
    W_E.set_row(2,{0,0,1});


    auto R1 = determinant(U_E*W_E* transpose(V_E)) * U_E*W_E* transpose(V_E);
    auto R2 = determinant(U_E* transpose(W_E)* transpose(V_E)) * U_E* transpose(W_E)* transpose(V_E);
    auto t1 = U_E.get_column(U_E.cols() - 1);
    auto t2 = -1* U_E.get_column(U_E.cols() - 1);

    auto point_candidates = point_options(K, points_0, points_1, R1, t1, R2, t2);

    // searching for the best combination of R and t
    std::vector<int> options_score;
    for (const auto &points :point_candidates){
        int option = 0;
        for (const auto &point: points){
            if (point.z()>0 ){
                option +=1;
            }
        }
        std::cout<<option<<std::endl;
        options_score.emplace_back(option);
    }


/*  ADDITION CODE for calculating in front of both cameras
    // searching for the best combination of R and t
    std::vector<int> options_score;
    for (const auto &points :point_candidates){  /// array of P vector3d of KPRT using point_0 and point_1
        int option = 0;
        for (const auto &point: points){
            if (point_0.z()>0 && point_1.z()>0){    // in front of two cameras
                option +=1;
            }
        }
        options_score.emplace_back(option);
    }
*/

    /*  CODE LEON
     std::cout << M0 << std::endl;
     int counter = 0;
     for (int i = 0; i < points_0.size(); i++){
         auto& p0 = points_0[i];
         auto& p1 = points_1[i];

         auto P = triangulate(p0, p1, M0, M1);

         Vector3D Q = R * P + t;
         if (P.z() > 0 && Q.z() > 0){
             counter++;
         }
         points.push_back(P);
     }
     std::cout << "In front: " << counter << std::endl;
     if (counter >= max_counter){
         max_counter = counter;
         best = pair;
         points_3d = points;
     }
 }
 return best;
  */


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

    //option 1 and 3 are preferred

    //TODO : VALIDATION;

    // we have a R and T and K for the M matrix
    // validation MP = p
    // pM = 0
    // projection matrix M

//    Matrix r_t3dcoor = Matrix(3,4);
//    r_t3dcoor.set_column(0,{R1.get_column(0)});  //set condition from above which R is the qualified R
//    r_t3dcoor.set_column(1,{R1.get_column(1)});  //set condition from above which R is the qualified R
//    r_t3dcoor.set_column(2,{R1.get_column(2)});  //set condition from above which R is the qualified R
//    r_t3dcoor.set_column(3,{t1});                          //set condition from above which t is the qualified t
//    auto Mproj = K*r_t3dcoor;
//
//    Vector pcomputed = mult(Mproj,points_3d);
//    //std::cout<<points_3d<<std::endl;
//
//    pcomputed - points_0



    //triangulate all corresponding image points (maybe we have to start with 8 instead of all the points)


/*
    std::vector<Vector3D> points_3d;
    // calculating 4d P point
    for (int i = 0; i < points_0.size(); i++) {

        auto x = points_0[i][0];
        auto y = points_0[i][1];
        auto x_prime = points_1[i][0];
        auto y_prime = points_1[i][1];

        Matrix identit = Matrix(3, 4);
        identit.set_row(0, {1, 0, 0, 0});
        identit.set_row(1, {0, 1, 0, 0});
        identit.set_row(2, {0, 0, 1, 0});

        auto M = K * identit;
        Matrix44 A;
        A.set_row(0, {x * M.get_row(2) - M.get_row(0)});
        A.set_row(1, {y * M.get_row(2) - M.get_row(1)});
        A.set_row(2, {x_prime * Mproj.get_row(2) - Mproj.get_row(0)});
        A.set_row(3, {y_prime * Mproj.get_row(2) - Mproj.get_row(1)});

        Matrix U = Matrix(A.rows(), A.cols(), 0.0);
        Matrix S = Matrix(A.rows(), A.cols(), 0.0);
        Matrix V = Matrix(A.cols(), A.cols(), 0.0);

        svd_decompose(A, U, S, V);
        //LEON Vector4D P4 = V.get_column(V.cols() - 1);
        //LEON return P4.cartesian();
        Vector3D p = V.get_column(V.cols() - 1);    //3d points calculated
        points_3d.push_back(p);
    }
    Vector points_2d = mult(Mproj,points_3d);     //validation of projected points are from Mproj*points_3d
    Vector nullaprox = points_2d * Mproj; //near to zero vector

*/

    return !points_3d.empty();
}
