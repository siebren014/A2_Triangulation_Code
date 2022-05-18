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



/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */
bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{

    if (!isvalid(points_0, points_1)){
        throw std::invalid_argument( "invalid input" );
    }

    // TODO: Estimate relative pose of two views. This can be subdivided into

    //Normalize points: origin should be at centre and average distance to centre should be sqrt2
    std::vector<Vector2D> points_0_normalized;
    std::vector<Vector2D> points_1_normalized;
    Vector2D sum_points0;
    Vector2D sum_points1;
    Vector2D sum_points0norm;
    Vector2D sum_points1norm;
    double dis0 = 0;
    double dis1 = 0;

    for (const auto& p1: points_0){
        sum_points0 += p1;
    }
    for (const auto& p1: points_1){
        sum_points1 += p1;
    }

    Vector2D mean_points0 = sum_points0/points_0.size();
    Vector2D mean_points1 = sum_points0/points_1.size();

    for (const auto& p1: points_0){
        dis0 += distance(p1, mean_points0);
    }
    auto avg_dis1 = dis0/points_0.size();
    for (const auto& p1: points_1){
        dis1 += distance(p1, mean_points1);
    }
    auto avg_dis2 = dis1/points_1.size();

    // if the average distance is bigger than square root of 2,normalize the points
//    if (avg_dis1 > sqrt(2)){
        auto norm_factor = avg_dis1/ sqrt(2);
        std::cout<<"avg_dis1 = "<<avg_dis1<< "norm_factor = "<< norm_factor<<std::endl;
        for (auto& p1: points_0){
            points_0_normalized.emplace_back(p1/norm_factor);
            sum_points0norm += p1/norm_factor;

        }
//    }
//    if (avg_dis2 > sqrt(2)){
        auto norm1_factor = avg_dis2/ sqrt(2);
        for (auto& p1: points_1){
            std::cout<<"p1 = "<<p1<< "p1/norm1_factor =  "<< p1/norm1_factor<<std::endl;
            points_1_normalized.emplace_back(p1/norm1_factor);
            sum_points1norm += p1/norm1_factor;
        }
//    }

    Vector2D mean_norm_points0 = sum_points0norm/points_0.size();
    Vector2D mean_norm_points1 = sum_points1norm/points_1.size();


    for (const auto& p1: points_0_normalized){
        dis0 += distance(p1, mean_norm_points0);
    }
    auto avg_dis0norm = dis0/points_0.size();
    for (const auto& p1: points_1_normalized){
        dis1 += distance(p1, mean_norm_points1);
    }
    auto avg_dis1norm = dis1/points_1.size();

    std::cout<<"dis 0= "<<avg_dis0norm/points_1.size() <<std::endl; // TODO: this should be sqrt 2 but oddly enough it isn't
    std::cout<<"dis 1= "<<avg_dis1norm/ points_1.size()<<std::endl;// TODO: this should be sqrt 2 but oddly enough it isn't


    /// define W_matrix based on amount of inputpoints
    Matrix W_matrix(points_0.size(), 9, 0.0);
    Matrix W_matrix_homo(points_0.size(), 9, 0.0);
    ///fill W_matrix by traversing through all the points.
    for (int i = 0; i < points_0.size(); i++){
        Vector2D p1 = points_0_normalized[i]; p1.homogeneous();
        Vector2D p2 = points_1_normalized[i]; p2.homogeneous();
        auto u1 = p1[0]; auto v1 = p1[1];
        auto u2 = p2[0]; auto v2 = p2[1];

        W_matrix.set_row(i, {points_0[i][0]*points_1[i][0], points_0[i][1]*points_1[i][0], points_1[i][0], points_0[i][0]*points_1[i][1], points_0[i][1]*points_1[i][1], points_1[i][1], points_0[i][0], points_0[i][1], 1});
        W_matrix_homo.set_row(i, {u1*u2, v1*u2, u2, u1*v2, v1*v2, v2, u1, v1, 1});
    }
    ///print W_matrix
    std::cout<<W_matrix<<std::endl;

    ///print W_matrix_homo
    std::cout<<W_matrix_homo<<std::endl;

    //use svd decompose to construct fundamental matrix from W matrix

    int m = R.rows()*2;
    int n = 12;

    Vector F = Vector(n, 0.0);
    Matrix U = Matrix(m,m,0.0);
    Matrix S = Matrix(m,n, 0.0);
    Matrix V = Matrix(n,n,0.0);

    svd_decompose(W_matrix, U, S, V);
    for (int i = 0; i < n; i++) {
        F[i] = V[i][n-1];

        std::cout<<F<<std::endl;

    //constraint matrix E


    // TODO:       - compute the essential matrix E;

    //    essential matrix = E = [T×]R

    ///     TODO:     - recover rotation R and t.
    //We can recover the R and t matrix from the fundamental matrix.
//    encodes information
//    about the camera matrices K,K′ and the relative translation T and rotation R between
//    the cameras.

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

    return points_3d.size() > 0;
}