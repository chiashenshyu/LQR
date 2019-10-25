#include <bits/stdc++.h>
#include "matplotlibcpp.h"
#include "d_are.hpp"
#pragma once
#include <cmath>
#include <string>
#define PLOT (1)

 
using namespace std;
//using namespace Eigen;
namespace plt = matplotlibcpp;

inline double mod2pi(double theta){
    return theta - 2*M_PI*floor(theta/2/M_PI); 
}

inline double pi2pi(double theta){
    return fmod(theta + M_PI, 2 * M_PI) - M_PI; 
}


struct state
{
    double x;
    double y;
    double w;
    double v;
    double L = 0.5;

    double steer_max = 45.0 / 180.0 * M_PI; 

    state(){x = 0;y = 0;w = 0;v = 0;}
    state(double x_,double y_, double w_,double v_){
        x = x_; y = y_; w = w_; v = v_;
    }
    //operator=(state& veh_)
    double distance (double xd , double yd){
        double dx = x-xd;
        double dy = y-yd;
        return sqrt(pow(dx,2)+pow(dy,2));
    }
    void update(double a,double delta,double dt){

    delta = (delta >  steer_max)?  steer_max : delta; 
    delta = (delta < -steer_max)? -steer_max : delta;   

    x =  x +  v * cos(w) * dt;
    y =  y +  v * sin(w) * dt;
    w =   w +  v / L * tan(delta) * dt;
    w = mod2pi(w);
    v =  v + a * dt;
    }
};


class LQRSteer {
    private:
    
    int pindx = -1;
    double Lfc = 2/*, k =0.8*/,L = 0.5, kv = 1,dt = 0.1;
    double target_speed = 10.0 / 3.6; // [m/s]
    double e, pe, th_e, pth_e;
    bool isRight = false; 

    vector<double> cy;
    vector<double> cx;
    vector<double> cd;
    vector<double> cdd;
    public:

    state veh;
    LQRSteer(vector<double> cx_,vector<double> cy_,vector<double> cd_,vector<double> cdd_,state veh_);
    vector<double> lqr_control(vector<double> sp,MatrixXd Q,MatrixXd R);
    double pid_vel();
    int calc_lookahead_pt();
    int last_target_idx = 0;

};
