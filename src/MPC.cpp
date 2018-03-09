#include "MPC.h"

using CppAD::AD;
using std::vector;
using namespace std;


void FG_polynomial::operator()(ADvector& fg, const ADvector& x) {
    fg[0] = 0.0;
    for (unsigned int j = 0; j < ptsx_.size(); j++) {
        // Observed point
        AD<double> y = 0;
        // Evaluate polynomial for ptsx(j)
        for (unsigned int i = 0; i < x.size(); i++) {
            y = y + x[i] * CppAD::pow(ptsx_[j], i);
        }
        fg[0] = fg[0] + CppAD::pow((ptsy_[j] - y), 2);
    }
}

// Implement global kinematic model
//x x_i=x_(i-1)+  cos⁡〖〖(psi〗_(i-1))* v_(i-1)*dt〗
//y y_i=y_(i-1)+  sin⁡〖〖(psi〗_(i-1))* v_(i-1)*dt〗
//v v_i=v_(i-1)+ a_(i-1)*dt
//psi 〖psi〗_i=〖psi〗_(i-1)+v_(i-1)  *  d_(i-1)*dt/Lf
//cte 〖cte〗_i=〖( polyeval(〗⁡〖x_(i-1))- y_(i-1)  )+〗  sin⁡〖〖(epsi〗_(i-1))* v_(i-1)*dt〗
//epsi  〖epsi〗_i=〖( 〖psi〗_(i-1 )- atan⁡polyeval'⁡〖x_(i-1) 〗 〗⁡〖 )+〗 v_(i-1)  *  d_(i-1)*dt/Lf

ADvector FG_motion::globalKinematic(const ADvector& state, const ADvector& actuators, const double& dt) {
    ADvector next_state(6);

    next_state[IDX_x] = state[IDX_x] + state[IDX_v] * CppAD::cos(state[IDX_psi]) * dt;
    next_state[IDX_y] = state[IDX_y] + state[IDX_v] * CppAD::sin(state[IDX_psi]) * dt;
    next_state[IDX_v] = state[IDX_v] + actuators[IDX_a] * dt;
    next_state[IDX_psi] = state[IDX_psi] - state[IDX_v] * actuators[IDX_delta] * dt / Lf;

    // Predicted y coordinate and angle
    AD<double> py = 0;
    AD<double> ppsi = 0;
    // Evaluate polynomial at x
    for (unsigned int i = 0; i < coeffs_.size(); i++) {
        py += coeffs_[i] * CppAD::pow(state[IDX_x], i);
    }
    // Evaluate Polynomial prime at x
    for (unsigned int i = 1; i < coeffs_.size(); i++) {
        ppsi += coeffs_[i] * CppAD::pow(state[IDX_x], i - 1) * i;
    }

    next_state[IDX_cte] = (py - state[IDX_y]) + (state[IDX_v] * CppAD::sin(state[IDX_epsi]) * dt);
    next_state[IDX_epsi] = (state[IDX_psi] - ppsi) - state[IDX_v] * actuators[IDX_delta] * dt / Lf;

    return next_state;
}

void FG_motion::operator()(ADvector& fg, const ADvector& x) {
    // Step 1: calculate Kinematic model
    // initial state is in attribute state
    // (to be optimized) actuators are in x
    // x is set of deltas and set of acceleration

    if (DEBUG) {
        cout << "__FG_motion()__" << endl;
        cout << "coefficients: " << coeffs_ << endl;
        cout << "Initial state: " <<state_ << endl;
        cout << "Initial actuators: " << actuators_ << endl;
        cout << "x: " << x << endl;
        cout << "fg: " << fg << endl;        
        cout << x.size() << endl;
        cout << fg.size() << endl;
    }

    const unsigned int time_steps = N;

    ADvector current_state(6);
    ADvector current_actuators(2);
    ADvector previous_actuators(2);
    ADvector next_state(6);

        // Cst to minimize
    fg[0] = 0.0;

    // Initialize current state to initial state given
    // Converting DVector to ADVector
    for (unsigned int i=0; i<state_.size(); i++) {
        current_state[i] = state_[i];
    }
    previous_actuators[IDX_delta] = actuators_[IDX_delta];
    previous_actuators[IDX_a]     = actuators_[IDX_a];


    for (unsigned int i = 0; i<time_steps; i++) {
         // Get delta and acceleration from optimizer
        current_actuators[IDX_delta] = x[i];
        current_actuators[IDX_a]     = x[i + time_steps];

        // determine next state based on current state and kinematic model
        // next_state = globalKinematic(current_state, current_actuators, DT);
        // relationship of current state + actuations and next state
        // based on our kinematic model
        next_state = globalKinematic(current_state, current_actuators, DT);

        // Update Cost
        // Penalize cross track error
        fg[0] += w_cte     * CppAD::pow(next_state[IDX_cte] - ref_cte, 2);
        // Penalize orientation
        fg[0] += w_epsi    * CppAD::pow(next_state[IDX_epsi] - ref_epsi, 2);
        // Penalize slow velocity
        fg[0] += w_v       * CppAD::pow(next_state[IDX_v] - ref_v, 2);


        // Penalize too much steering and acceleration
        fg[0] += w_delta   * CppAD::pow(current_actuators[IDX_delta], 2);
        fg[0] += w_a       * CppAD::pow(current_actuators[IDX_a], 2);

        // Penalize change in steering and acceleration
        fg[0] += w_delta_t * CppAD::pow(previous_actuators[IDX_delta] - current_actuators[IDX_delta], 2);
        fg[0] += w_a_t     * CppAD::pow(previous_actuators[IDX_a] - current_actuators[IDX_a], 2);

        // Set x coordinates
        fg[1 + i] = next_state[IDX_x] - x[2*time_steps + i];
        fg[1 + i + time_steps] = next_state[IDX_y] - x[3*time_steps + i];

        // Repeat with new state
        current_state = next_state;
        previous_actuators = current_actuators;


    }

    if (DEBUG)
        cout << "Cost #" << fg[0] << endl;


}

MPC::MPC() {

    // Initialize coefficients
    coeffs.resize(DEGREE_POLYNOMIAL);
    coeffs_x_lowerbound.resize(DEGREE_POLYNOMIAL);
    coeffs_x_upperbound.resize(DEGREE_POLYNOMIAL);
    // No constraints for coefficients
    for (unsigned int i = 0; i < coeffs.size(); i++) {
        coeffs_x_lowerbound[i] = -1e19;
        coeffs_x_upperbound[i] = 1e19;
    }

    //**************************************************************
    //* SET INITIAL VALUES OF VARIABLES
    //**************************************************************
    x.resize(NX);
    x_lowerbound.resize(NX);
    x_upperbound.resize(NX);

    // Steering limits set to 25 radians
    for (int i = 0; i < N; i++) {
        x[i] = 0;
        x_lowerbound[i] = -0.436332;
        x_upperbound[i] = 0.436332;
    }
    // throttle limits
    for (int i = N; i < 2*N; i++) {
        x[i] = 0.1; // Initialize to slight acceleration / thottle
        x_lowerbound[i] = -1.0;
        x_upperbound[i] = 1.0;
    }    
    // x and y limits
    for (int i = 2*N; i < NX; i++) {
        x[i] = 0.;
        x_lowerbound[i] = -1e19;
        x_upperbound[i] = 1e19;
    }

    g_lowerbound.resize(2*N);
    g_upperbound.resize(2*N);
    for (int i = 0; i < 2*N; i++) {
        g_lowerbound[i] = 0;
        g_upperbound[i] = 0;
    }    

    next_xs.resize(N);
    next_ys.resize(N);
    mpc_xs.resize(N);
    mpc_ys.resize(N);

    // object that computes objective and constraints
    // Uncomment this if you'd like more print information
    if (DEBUG)
        options_model += "Integer print_level  5\n";
    else
        options_model += "Integer print_level  0\n";
    // NOTE: Setting sparse to true allows the solver to take advantage
    // of sparse routines, this makes the computation MUCH FASTER.
    options_model += "Sparse  true        forward\n";
    options_model += "Sparse  true        reverse\n";
    // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
    options_model += "Numeric max_cpu_time         30.0\n";    // object that computes objective and constraints
    // Uncomment this if you'd like more print information
    options_poly += "Integer print_level  0\n";
    // NOTE: Setting sparse to true allows the solver to take advantage
    // of sparse routines, this makes the computation MUCH FASTER.
    options_poly += "Sparse  true        forward\n";
    options_poly += "Sparse  true        reverse\n";
    // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
    options_poly += "Numeric max_cpu_time         30.0\n";
}

MPC::~MPC() {};

void MPC::solve(const Dvector& state, const Dvector& actuators) {
    if (DEBUG)
        cout << "__solve__ enter" << endl;

    if (LOG) {
        log.open("logfile.txt",std::ios_base::app | std::ios_base::out);
        log << "state:\n";
        for (unsigned int i = 0;i<state.size();i++)
            log << state[i] << "\n";
        log << "actuators:\n";
        for (unsigned int i = 0;i<actuators.size();i++)
            log << actuators[i] << "\n";
        log.close();
    }

    CppAD::ipopt::solve_result<Dvector> solution;
    bool ok = true;

    FG_motion fg_motion(coeffs, state, actuators);

   // solve the problem
    CppAD::ipopt::solve<Dvector, FG_motion>(
     options_model,
     x,
     x_lowerbound,
     x_upperbound,
     g_lowerbound,
     g_upperbound,
     fg_motion,
     solution);

    auto cost = solution.obj_value;

    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

    if (ok) {
        cout << "OK! Cost: " << cost << endl;
        if (DEBUG) {
            cout << "OK! Solution x: " << solution.x << endl;
            cout << "OK! Solution g: " << solution.g << endl;
        }
    }
    else {
        cout << "SOMETHING IS WRONG! " << solution.status << endl;
    }

    // To cater for latency, ignore first value at index 0
    steer = solution.x[0];
    throttle = solution.x[N + 0];

    mpc_xs.clear();
    mpc_ys.clear();
    for (unsigned int i =0; i<N;i++) {
        mpc_xs.push_back(solution.x[2*N + i]);
        mpc_ys.push_back(solution.x[3*N + i]);
    }

    if (DEBUG)
        cout << "Size of solution x: " << x.size() << endl;
}

void MPC::setPolynomialPoints() {

    next_xs.clear();
    next_ys.clear();

    for (unsigned int i = 0; i < N; ++i) {
        const double dx = D * i;
        next_xs.push_back(dx);
        next_ys.push_back(polyeval(dx));
    }
}


// // Adapted from lecture notes:
// // Estimate next state based on current state
Dvector MPC::globalKinematic(const Dvector& state, const Dvector& actuators, double dt) {
    // Create a new vector for the next state.
    Dvector next_state(6);

    next_state[IDX_x] = state[IDX_x] + state[IDX_v] * cos(state[IDX_psi]) * dt;
    next_state[IDX_y] = state[IDX_y] + state[IDX_v] * sin(state[IDX_psi]) * dt;
    next_state[IDX_v] = state[IDX_v] + actuators[IDX_a] * dt;
    next_state[IDX_psi] = state[IDX_psi] - state[IDX_v] * actuators[IDX_delta] * dt / Lf;

    // Predicted y coordinate and angle
    double py = 0;
    double ppsi = 0;
    // Evaluate polynomial at x
    for (unsigned int i = 0; i < coeffs.size(); i++) {
        py += coeffs[i] * pow(state[IDX_x], i);
    }
    // Evaluate Polynomial prime at x
    for (unsigned int i = 1; i < coeffs.size(); i++) {
        ppsi += coeffs[i] * pow(state[IDX_x], i - 1) * i;
    }

    next_state[IDX_cte] = (py - state[IDX_y]) + (state[IDX_v] * sin(state[IDX_epsi]) * dt);
    next_state[IDX_epsi] = (state[IDX_psi] - ppsi) - state[IDX_v] * actuators[IDX_delta] * dt / Lf;

    return next_state;
}

// Fit polynomial to Waypoints
// Reset to car's point of view

void MPC::fitWaypoints(vector<double>& ptsx, vector<double>& ptsy) {
    // Solve polynomial

    CppAD::ipopt::solve_result<Dvector> solution;
    bool ok = true;

    // Initialize values
    for (unsigned int i = 0; i < coeffs.size(); i++) {
        coeffs[i] = 1;
    }
    
    // Map vector values to Devcor values
    Dvector ptsx_d;
    Dvector ptsy_d;
    ptsx_d.resize(ptsx.size());
    ptsy_d.resize(ptsy.size());
    for (unsigned int i = 0; i<ptsx.size(); i++){
        ptsx_d[i] = ptsx[i];
        ptsy_d[i] = ptsy[i];
    }
    FG_polynomial fg_polynomial(ptsx_d, ptsy_d);

    CppAD::ipopt::solve<Dvector, FG_polynomial>(
        options_poly,
        coeffs,
        coeffs_x_lowerbound,
        coeffs_x_upperbound,
        coeffs_g_lowerbound,
        coeffs_g_upperbound,
        fg_polynomial,
        solution);

    auto cost = solution.obj_value;

    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;
    if (ok) {
        cout << "OK! Polynomial Cost: " << cost << endl;
        for (unsigned int i = 0; i < coeffs.size(); i++) {
            coeffs[i] = solution.x[i];
        }
        if (DEBUG)
            cout << coeffs << endl;

    }
    else {
        cout << "SOMETHING IS WRONG in Polynomial ! " << solution.status << endl;
    }

    setPolynomialPoints();
}

// Evaluate a polynomial.
double MPC::polyeval(double x) {
    double result = 0.0;
    for (unsigned int i = 0; i < coeffs.size(); i++) {
        result += coeffs[i] * pow(x, i);
    }
    return result;
}

double MPC::polyeval_t(double x) {
    double result = 0.0;
    for (unsigned int i = 1; i < coeffs.size(); i++) {
        result += coeffs[i] * pow(x, i - 1) * i;
    }
    return result;
}
