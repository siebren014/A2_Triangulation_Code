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

    Vector2D mean_points0 = sum_points0 / points_0.size();
    Vector2D mean_points1 = sum_points0 / points_1.size();
    for (const auto &p1: points_0) {
        dis0 += distance(p1, mean_points0);
    }
    auto avg_dis1 = dis0 / points_0.size();
    for (const auto &p1: points_1) {
        dis1 += distance(p1, mean_points1);
    }
    auto avg_dis2 = dis1 / points_1.size();

    auto norm_factor = avg_dis1 / sqrt(2);
    for (auto &p1: points_0) {
        points_0_normalized.emplace_back(p1 / norm_factor);
        sum_points0norm += p1 / norm_factor;

    }

    auto norm1_factor = avg_dis2 / sqrt(2);
    for (auto &p1: points_1) {
        points_1_normalized.emplace_back(p1 / norm1_factor);
        sum_points1norm += p1 / norm1_factor;
    }

    Vector2D mean_norm_points0 = sum_points0norm / points_0.size();
    Vector2D mean_norm_points1 = sum_points1norm / points_1.size();

    for (const auto &p1: points_0_normalized) {
        dis0_norm += distance(p1, mean_norm_points0);
    }
    for (const auto &p1: points_1_normalized) {
        dis1_norm += distance(p1, mean_norm_points1);
    }
    normalize_results.emplace_back(points_0_normalized);
    normalize_results.emplace_back(points_1_normalized);

    return normalize_results;
}
std::vector<Vector3D> Points(const Matrix33 &K, const std::vector<Vector2D> &points_0,const std::vector<Vector2D> &points_1, const Matrix &R, const Vector &t){

    Matrix k_t = Matrix(3,4);
    k_t.set_column(0,{R.get_column(0)});
    k_t.set_column(1,{R.get_column(1)});
    k_t.set_column(2,{R.get_column(2)});
    k_t.set_column(3,{t});


    Matrix identit = identity(4,4);

    auto M_prime = K*k_t;
    auto M = K* identit;

    std::vector<Vector3D> points_3d;

    for (int i = 0; i < points_0.size(); i++){

        auto x = points_0[i][0];
        auto y = points_0[i][1];
        auto x_prime = points_1[i][0];
        auto y_prime = points_1[i][1];

        Matrix A;
        A.set_row(0,{x*M.get_row(3) - M.get_row(1)});
        A.set_row(1,{y*M.get_row(3) - M.get_row(2)});
        A.set_row(2,{x_prime*M.get_row(3) - M_prime.get_row(1)});
        A.set_row(3,{y_prime*M.get_row(3) - M_prime.get_row(1)});

        std::cout<<"A"<<A<<std::endl;
    }

    return points_3d;
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


    /// define W_matrix based on amount of inputpoints
    Matrix W_matrix(points_0.size(), 9, 0.0);
    Matrix W_matrix_homo(points_0.size(), 9, 0.0);
    ///fill W_matrix by traversing through all the points.
    for (int i = 0; i < points_0.size(); i++){
        Vector2D p1 = points_0_norm[i];
        Vector2D p2 = points_1_norm[i];
        auto u1 = p1[0]; auto v1 = p1[1];
        auto u2 = p2[0]; auto v2 = p2[1];

        W_matrix.set_row(i, {points_0[i][0]*points_1[i][0], points_0[i][1]*points_1[i][0], points_1[i][0], points_0[i][0]*points_1[i][1], points_0[i][1]*points_1[i][1], points_1[i][1], points_0[i][0], points_0[i][1], 1});
        W_matrix_homo.set_row(i, {u1*u2, v1*u2, u2, u1*v2, v1*v2, v2, u1, v1, 1});
    }

    int m = points_0.size();
    int n = 9;

//    Vector F = Vector(n, 0.0);
    Matrix U = Matrix(m,m,0.0);
    Matrix S = Matrix(m,n, 0.0);
    Matrix V = Matrix(n,n,0.0);

    svd_decompose(Matrix (W_matrix_homo), U, S, V);
    Vector F = V.get_column(V.cols() - 1);


    std::cout<<F<<std::endl;
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

    //set the rank to 2
    S_mat[2][2]=0;


    Matrix33 F_bestrank = (U_mat * S_mat * V.transpose());


    /// computed assignment. Divide each element by v

    auto F_scaled = F_bestrank/ F_bestrank[2][2];


    // TODO: STEP 2

    Matrix33 K;
    K.set_row(0,{fx, 0, cx});
    K.set_row(1,{0,fy,cy});
    K.set_row(2,{0,0,1});

    std::cout<<"K"<<K<<std::endl;

    //E Matrix
    Matrix E = transpose(K)*F_mat*K;

    Matrix U_E = Matrix(E.rows(),E.rows(),0.0);
    Matrix V_E = Matrix(E.rows(),E.cols(),0.0);
    Matrix S_E = Matrix(E.cols(),E.cols(),0.0);

    svd_decompose(Matrix (E), U_E, S_E, V_E);

    // R and t have 2 potential values, so 4 values. Means that we have 4 candidates.
    Matrix W_E = Matrix(3,3);
    W_E.set_row(0,{0,-1,0});
    W_E.set_row(1,{1,0,0});
    W_E.set_row(2,{0,0,1});

    std::cout<<"matrix U = "<<U_E<<std::endl;
    std::cout<<"matrix v = "<<V_E<<std::endl;
    std::cout<<"matrix w = "<<W_E<<std::endl;

    auto R1 = determinant(U_E*W_E* transpose(V_E)) * U_E*W_E* transpose(V_E);
    auto R2 = determinant(U_E* transpose(W_E)* transpose(V_E)) * U_E* transpose(W_E)* transpose(V_E);
    auto t1 = U_E.get_column(U_E.cols() - 1);
    auto t2 = -1* U_E.get_column(U_E.cols() - 1);

    auto option_1 = Points(K, points_0, points_1, R1, t1);
    auto option_2 = Points(K, points_0, points_1, R1, t2);
    auto option_3 = Points(K, points_0, points_1, R2, t1);
    auto option_4 = Points(K, points_0, points_1, R2, t2);


//
//
//    Matrix k_t = Matrix(3,4);
//    k_t.set_column(0,{R1.get_column(0)});
//    k_t.set_column(1,{R1.get_column(1)});
//    k_t.set_column(2,{R1.get_column(2)});
//    k_t.set_column(3,{t1});
//
//
//    Matrix identit = identity(4,4);
//
//    auto M_prime = K*k_t;
//    auto M = K* identit;
//
//    auto x =(M.get_row(3)*P) - (M.get_row(1)*P);
//    auto y =(M.get_row(3)*P) - (M.get_row(2)*P);
//    auto z =(M.get_row(2)*P) - (M.get_row(1)*P);
//
//    Matrix A;
//    A.set_row(0,{x*M.get_row(3) - M.get_row(1)})
//    A.set_row(1,{y*M.get_row(3) - M.get_row(2)})
//    A.set_row(2,{x_prime*M.get_row(3) - M_prime.get_row(1)})
//    A.set_row(3,{y_prime*M.get_row(3) - M_prime.get_row(1)})

//
//    std::cout<<t1<<std::endl;


//    std::cout<<r_test<<std::endl;
//    auto M = K*(R1*t1);
////    auto x(M3*P) - (M1*P) = 0;
//





    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       There are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.

    return !points_3d.empty();
}
