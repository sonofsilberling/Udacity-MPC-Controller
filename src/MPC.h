#ifndef MPC_H
#define MPC_H

#include <cppad/ipopt/solve.hpp>
#include <vector>
#include <ios>
#include <fstream>

using namespace std;
using CppAD::AD;

const bool DEBUG = false;
const bool LOG = false;

// Time between predictions ahead
const double DT = 0.1; 
const double Lf = 2.67;

const int N = 10; // how many states we "lookahead" into the future

const int NUMBER_OF_STATES = 6; // px, py, psi, v, cte, epsi
const int NUMBER_OF_ACTUATORS = 2; // steering angle, acceleration

const int DEGREE_POLYNOMIAL = 4;

// number of state + actuation variables
// EXCLUDING initial state
// How many variables the optimizer has to concern itself with (x)
// N steering values
// N throttle values. AND
// N x values - required for output
// N y values - required for output but not strictly the optimizer
const int NX = N * 4 ;
// All input values for the optimizer
const int NG = N * ( NUMBER_OF_STATES );

// Index in State vector
const int IDX_x = 0;
const int IDX_y = 1;
const int IDX_psi = 2;
const int IDX_v = 3;
const int IDX_cte = 4;
const int IDX_epsi = 5;
// Index in actuator vector
const int IDX_delta = 0;
const int IDX_a = 1;

// Reference CTE, PSI, and target velocity
const double ref_cte = 0;
const double ref_epsi = 0;
const double ref_v = 100;

// Weights for cost function
const double w_cte = 2000;
const double w_epsi = 2000;
const double w_v = 1;
const double w_delta = 5;
const double w_a = 1;
const double w_delta_t = 200;
const double w_a_t = 10;

typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
typedef CPPAD_TESTVECTOR(double) Dvector;

// IPopt functor for motion model
class FG_motion {
public:
    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
    // Coefficients of fitted polynomial curve
    Dvector coeffs_;

    // Initial (constant) state
    Dvector state_;
    Dvector actuators_;

    FG_motion(Dvector coeffs_in, Dvector state_in, Dvector actuators_in) : coeffs_(coeffs_in), state_(state_in), actuators_(actuators_in) {
    }

    void operator()(ADvector& fg, const ADvector& x);

    ADvector globalKinematic(const ADvector& state, const ADvector& actuators, const double& dt);
};

// IPopt functor for polynomial /curve fitting
class FG_polynomial {
public:
    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
    Dvector ptsx_;
    Dvector ptsy_;

    FG_polynomial(Dvector ptsx_in, Dvector ptsy_in) : ptsx_(ptsx_in), ptsy_(ptsy_in) {
    }

    void operator()(ADvector& fg, const ADvector& x);
};

class MPC {
public:

    // Distance for predicted waypoints
    const double D = 5.0;

    std::ofstream log;

    // Predicated Steer and Throttle values
    double steer;
    double throttle;

    Dvector x; // where all the state and actuation variables will be stored
    Dvector x_lowerbound; //lower limit for each corresponding variable in x
    Dvector x_upperbound; //upper limit for each corresponding variable in x
    Dvector g_lowerbound; // value constraint for each corresponding constraint expression
    Dvector g_upperbound; // value constraint for each corresponding constraint expression

    // Predicted waypoints based on polynomial - always N of those
    vector<double> next_xs;
    vector<double> next_ys;
    // Predicted waypoints based on optimizer - always N of those
    vector<double> mpc_xs;
    vector<double> mpc_ys;

    // Fitted polynomial coefficients
    Dvector coeffs;
    Dvector coeffs_x_lowerbound;
    Dvector coeffs_x_upperbound;
    Dvector coeffs_g_lowerbound;
    Dvector coeffs_g_upperbound;

    // options for IPOPT solver
    std::string options_poly;
    std::string options_model;

    MPC();
    virtual ~MPC();

    Dvector globalKinematic(const Dvector& state, const Dvector& actuators, double dt);

    void setPolynomialPoints();
    void fitWaypoints(vector<double>& ptsx, vector<double>& ptsy);

    // Evaluate Polynomial
    double polyeval(double x);

    // Evaluate first derivative of polynomial
    double polyeval_t(double x);
    //
    // this function solves the model given the current state and road curve coefficients.
    void solve(const Dvector& state, const Dvector& actuators);

};

#endif /* MPC_H */
