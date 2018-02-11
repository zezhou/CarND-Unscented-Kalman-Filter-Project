#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd residual_sum = VectorXd(4);
  residual_sum  << 0, 0, 0, 0;
  for ( int i = 0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    residual_sum += residual;
  }
  VectorXd rmse = VectorXd(4);
  rmse << 0, 0, 0, 0;
  rmse = residual_sum / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}