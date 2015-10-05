#include <vector>
#include <algorithm> // max_element, max
#include <memory>
#include <limits>
#include "ODD_delayed_values.h"
//#include "ODD_rhs.h"
#include "ODD_stepper.h"
#include "ODD_utils.h" 

using namespace std;

ODD_stepper::ODD_stepper() {
}

CashKarpStepper::CashKarpStepper() {
}

void CashKarpStepper::init(ODD_rhs user_rhs, const ODD_parameter* pars, const ODD_solution* sol, const ODD_delays* tau, const vector<double>* p) {
    this->user_rhs = user_rhs;

    size_t Ndim = pars->find("N")->second;
    k1 = vector<double>(Ndim);
    k2 = vector<double>(Ndim);
    k3 = vector<double>(Ndim);
    k4 = vector<double>(Ndim);
    k5 = vector<double>(Ndim);
    k6 = vector<double>(Ndim);
    f0 = vector<double>(Ndim);
    x0 = vector<double>(Ndim);
    last_point = vector<double>(Ndim);
    x_tmp = vector<double>(Ndim);

    // minimal and maximal stepwidth
    hmin = pars->find("hmin")->second;
    hmax = pars->find("hmax")->second;

    h_next = min(pars->find("h0")->second, hmax);
    // relative and absolute tolerances
    rel_tol = pars->find("rel_tol")->second;
    abs_tol = pars->find("abs_tol")->second;
    // bounding stepsize adaptation scalefactors
    minimal_stepsize_rescale = pars->find("minimal_stepsize_rescale")->second;
    maximal_stepsize_rescale = pars->find("maximal_stepsize_rescale")->second;

    safety_factor = pars->find("safety_factor")->second;

    tau0_offset = pars->find("tau0_offset")->second;

    this->sol = sol;
    this->tau = tau;
    this->p = p;
}

void CashKarpStepper::setBasePoint(double t, const vector<double>& x, const vector<double>& f) {
    t0 = t;
    x0 = x;
    f0 = f;
}

bool CashKarpStepper::step() {
    return step(h_next);
}

bool CashKarpStepper::step(double h) {
    //    cout << "called step(h) with h = " << h << endl;
    h = min(h, hmax);

    bool success; // to be returned

    // remember the value of h
    h_last = h;

    /* calculate runge kutta coefficients */

    // calc k1
    k1 = h*f0;
    //    cout << "k1 = " << k1 << endl;

    // get delayed values for next runge kutta point
    sol->interpolate(t0 + h * bt.a2 - *tau, *xd);
    x_tmp = x0 + bt.b21*k1;
    xd->fillZeroDelays(x_tmp);
    // calc k2
    rhs(t0 + h * bt.a2, x_tmp, *xd, *p, k2);
    k2 = h*k2;

    //    cout<<"k2 = "<<k2<<endl;

    // get delayed values for next runge kutta point
    sol->interpolate(t0 + h * bt.a3 - *tau, *xd);
    x_tmp = x0 + bt.b31 * k1 + bt.b32*k2;
    xd->fillZeroDelays(x_tmp);
    // calc k3
    rhs(t0 + h * bt.a3, x_tmp, *xd, *p, k3);
    k3 = h*k3;

    //    cout<<"k3 = "<<k3<<endl;

    // get delayed values for next runge kutta point
    sol->interpolate(t0 + h * bt.a4 - *tau, *xd);
    x_tmp = x0 + bt.b41 * k1 + bt.b42 * k2 + bt.b43*k3;
    xd->fillZeroDelays(x_tmp);
    // calc k4
    rhs(t0 + h * bt.a4, x_tmp, *xd, *p, k4);
    k4 = h*k4;

    //    cout<<"k4 = "<<k4<<endl;

    // get delayed values for next runge kutta point
    sol->interpolate(t0 + h * bt.a5 - *tau, *xd);
    x_tmp = x0 + bt.b51 * k1 + bt.b52 * k2 + bt.b53 * k3 + bt.b54*k4;
    xd->fillZeroDelays(x_tmp);
    // calc k5
    rhs(t0 + h * bt.a5, x_tmp, *xd, *p, k5);
    k5 = h*k5;

    //    cout<<"k5 = "<<k5<<endl;

    // get delayed values for next runge kutta point
    sol->interpolate(t0 + h * bt.a6 - *tau, *xd);
    x_tmp = x0 + bt.b61 * k1 + bt.b62 * k2 + bt.b63 * k3 + bt.b64 * k4 + bt.b65*k5;
    xd->fillZeroDelays(x_tmp);
    // calc k6
    rhs(t0 + h * bt.a6, x_tmp, *xd, *p, k6);
    k6 = h*k6;

    //    cout<<"k6 = "<<k6<<endl;

    /* calculate solution estimates at t0 + h */

    // order 4 increment (use storage k2 because it isn't needed anymore)
    increment_order_4 = bt.c41 * k1 + bt.c43 * k3 + bt.c44 * k4 + bt.c46 * k6;
    // order 5 increment
    increment_order_5 = bt.c51 * k1 + bt.c53 * k3 + bt.c54 * k4 + bt.c55 * k5 + bt.c56 * k6;

    /* approximated solution value 
     * NOTE: uses order 5 instead order 4 as recommended by Puls and Heitsch.
     * (Shampine ("global order estimates for odes" or similar does this, too)) 
     * They say that using the order 4 solution is better for stiff problems.
     */
    last_point = x0 + increment_order_5;

    // componentwise error estimation
    abs_err = cwise_abs(increment_order_4 - increment_order_5) * pow(h, 4);
    double max_abs_err = *max_element(abs_err.begin(), abs_err.end());
    double max_rel_err;
    
    if (rel_tol > 0) {
        component_scales = cwise_abs(k1) + cwise_abs(x0) + EPS; // XXX: why multiplicate f0 with h here (k1==h*f0)? :leo: found (and lost) source with recommendation ("fractional errors stay constant"), but no explanation
        rel_err = cwise_divide(abs_err, component_scales);
        max_rel_err = *max_element(rel_err.begin(), rel_err.end());
        /* check whether this step fulfills tolerance requirements */
        success = max_rel_err < rel_tol && max_abs_err < abs_tol;
        /* adapt stepsize
         * XXX: that is NOT the same formula as Dimitry was using :leo: 
         */
        h_next = h * safety_factor * pow(min(rel_tol / max(max_rel_err, EPS), abs_tol / max(max_abs_err, EPS)), 0.25);
    } else {
        /* check whether this step fulfills tolerance requirements */
        success = max_abs_err < abs_tol;
        h_next = h * safety_factor * pow(abs_tol / max(max_abs_err, EPS), 0.25);
    }

    h_next = max(max(h_next, hmin), h * minimal_stepsize_rescale);
    h_next = min(min(h_next, hmax), h * maximal_stepsize_rescale);
    h_last = h;


    if (std::isnan(h_next)) {
        throw runtime_error("h_next is nan.");
    }


//    //    if (t0 > 10080) {
//    if (t0 > 5080 && t0 < 5100) {
//        cout << "t0 = " << t0 << endl;
//        auto prec = cout.precision();
//        cout.precision(numeric_limits<double>::digits10);
//        cout << "performed step with stepsize h = " << h << endl;
//        cout << "max_rel_err = " << max_rel_err << endl;
//        cout << "max_abs_err = " << max_abs_err << endl;
//        cout.precision(prec);
//
//        cout << "success: " << success << endl;
//    }
    return success;
}
