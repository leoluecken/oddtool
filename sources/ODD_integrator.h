/* 
 * File:   ODD_integrator.h
 * Author: luecken
 *
 * Created on 27. Januar 2014, 14:13
 * 
 * 
 */

#ifndef ODD_INTEGRATOR_H
#define	ODD_INTEGRATOR_H

#include "ODD_solution.h"
#include "ODD_stepper.h"
#include "ODD_output.h"
#include "ODD_input.h"
#include "ODD_utils.h"
#include "ODD_delayed_values.h"
#include <memory>
#include <string>
#include <cmath>
#include <map>
#include <vector>

/** ODD_integrator
 * 
 * The template class ODD_integrator<STEPPER> contains an abstract interface to
 * a DDE integrator with a custom stepper (should extend ODD_stepper), which performs
 * the single integration steps, including error estimation / stepsize control
 * 
 */
template<typename STEPPER>
class ODD_integrator {
    // name of the integrator (corresponds to the root of the filenames)
    std::string name;
    // system dimension
    int N;
    // minimal_delay
    double tau_min;
    // integration start time
    double t_begin;
    // integration end time
    double t_end;
    // maximal number of steps during the integration
    double max_steps;
    // nr of steps calculated by the stepper that did satisfy error bounds 
    int steps_success = 0;
    // nr of steps calculated by the stepper that did not satisfy error bounds 
    int steps_fail = 0;

    // user supplied history function
    ODD_history_fct* history_fct;
    // user supplied rhs function
    ODD_rhs* user_rhs;
    // wrapped rhs function

    void rhs(double t, const std::vector<double>& x, ODD_delayed_values& xd, const std::vector<double>& p, std::vector<double>& result) {
        // switch to user indices while calling user supplied rhs
        xd.useUserIndices(true);
        user_rhs(t, x, xd, p, result);
        xd.useUserIndices(false);
    }

    // integration parameters
    //    std::shared_ptr<const ODD_parameter> pars;
    const ODD_parameter* pars;
    // parameters for rhs function
    //    std::shared_ptr<const std::vector<double> > rhs_pars;
    const std::vector<double>* rhs_pars;
    // delay values
    //    std::shared_ptr<const ODD_delays> tau;
    const ODD_delays* tau;

    // input
    std::unique_ptr<ODD_input> input;
    // solution
    std::unique_ptr<ODD_solution> sol;
    // numeric stepper
    std::unique_ptr<STEPPER> stepper;

    // actual time
    double t;
    // actual state
    std::vector<double> x;
    // actual derivative
    std::vector<double> f;

    // temp storage for handling delayed values x(t-tau) of the solution 
    std::unique_ptr<ODD_delayed_values> xd;

    // check out the user supplied rhs function
    void inspect_rhs();

    // accept the last step taken by the stepper
    bool takeStep();

    // accept the last step taken by the stepper
    void acceptStep();

public:
    // constructor
    ODD_integrator(ODD_rhs rhs, ODD_history_fct hist = nullptr, std::string parameter_filename = "ODD_parameters.txt");
    void integrate();
    void integrate(double t_end_new);
};


// constructor initializes the members

template<typename STEPPER>
ODD_integrator<STEPPER>::ODD_integrator(ODD_rhs user_rhs, ODD_history_fct hist, std::string parameter_filename) {
    // create an input
    //    input = std::make_shared<ODD_input > ();
    input = std::unique_ptr<ODD_input > (new ODD_input());
    // load the data from the parameter file and history file(s) (specified in parameter file) into the input
    input->init(parameter_filename, hist);

    std::cout << "\ninput done.\n" << std::endl;

    // get parameters from input
    pars = input->getParameters();
    rhs_pars = input->getRHSParameters();
    tau = input->getDelays();
    name = input->getName();
    tau_min = pars->find("tau_min")->second;
    N = pars->find("N")->second;
    t_begin = pars->find("t_start")->second;
    t_end = pars->find("t_end")->second;
    max_steps = pars->find("max_steps")->second;
    history_fct = hist;
    this->user_rhs = user_rhs;

    // allocate x and f
    x = std::vector<double>(N);
    f = std::vector<double>(N);

    // init solution object
    //    sol = std::make_shared<ODD_solution > ();
    sol = std::unique_ptr<ODD_solution > (new ODD_solution());
    if (hist == nullptr || pars->find("force_history_from_file")->second != 0) {
        // if no history function supplied, data must have been provided in file
        sol->init(name, pars, input->getHistory());
    } else {
        // history function was given "analytically"
        sol->init(name, pars, hist);
    }

    // init data structure for passing around delayed values
    //    xd = std::make_shared<ODD_delayed_values > ();
    xd = std::unique_ptr<ODD_delayed_values > (new ODD_delayed_values());
    xd->init(input->getUserDelayIndices(), pars->find("tau0_offset")->second);
    inspect_rhs();

    std::cout << "\nstorage initialized. size = " << sol->size() << std::endl;

    // init the stepper
    //    stepper = std::make_shared<STEPPER > ();
    stepper = std::unique_ptr<STEPPER > (new STEPPER());
    stepper->init(user_rhs, pars, sol.get(), tau, rhs_pars);
    stepper->init_xd(xd.get()); // uses the same temp container
}





// integrate() generates a solution (this->sol) in the
// interval [sol->lastT(), t_end] using the stepper (this->stepper)

template<typename STEPPER>
void ODD_integrator<STEPPER>::integrate() {

//    double t_trace = 0;
//    size_t last_steps_success = 0;
//    size_t last_steps_fail = 0;





    std::cout << "\n\nODD_integrator::integrate() called." << std::endl;
    std::cout << "t_begin = " << t_begin << std::endl;
    std::cout << "t_end = " << t_end << std::endl;

    // get initial data
    t_begin = sol->lastT();
    t = t_begin;
    x = sol->lastX();

    // trace initial data
    std::cout << "t0 = " << t_begin << std::endl;
    std::cout << "x0 = [ ";
    for (double xi : x) std::cout << xi << " ";
    std::cout << "]" << std::endl;

    // trace tau
    std::cout << "tau = [ ";
    for (double taui : *tau) std::cout << taui << " ";
    std::cout << "]" << std::endl;

    // calculate initial derivative
    sol->interpolate(t - *tau, *xd);
    xd->fillZeroDelays(x);
    rhs(t, x, *xd, *rhs_pars, f);

    // assure that at least REQUIRED_INTERPOLATION_POINTS-1 steps are taken in [t_start, t_start + tau_min],
    // so no interpolation attempt is made on t[wi] when it still doesn't contain REQUIRED_INTERPOLATION_POINTS 
    // points (after initialization t[wi].size() == 1)
    if (tau_min != 0) {
        double orig_hmax = stepper->getHmax();
        double max_first_step = tau_min / (sol->requiredInterpolationPoints());
        stepper->setHmax(max_first_step);
        for (size_t k = 0; k != sol->requiredInterpolationPoints() - 1 && t < t_end; ++k) {
            bool success = takeStep();
            if (!success) {
                stepper->setHmax(orig_hmax);
                goto integration_end;
            }
            if (steps_success == max_steps) {
                stepper->setHmax(orig_hmax);
                std::clog << "WARNING: maximal number of steps reached.\nStopping integration at time t = " << t << std::endl;
                goto integration_end;
            }
        }
        stepper->setHmax(orig_hmax);
    }

    while (t < t_end) {
        //                std::cout << "t = " << t << std::endl;
        bool success = takeStep();


//        if (t > t_trace) {
//            std::cout << "steps_success = " << steps_success - last_steps_success << std::endl;
//            std::cout << "steps_fail = " << steps_fail - last_steps_fail << std::endl;
//            t_trace = t + 10000;
//            last_steps_success = steps_success;
//            last_steps_fail = steps_fail;
//        }

        if (!success) goto integration_end;
        if (steps_success == max_steps) {
            std::clog << "WARNING: maximal number of steps reached.\nStopping integration at time t = " << t << std::endl;
            break;
        }
    }

integration_end:
    //    std::cout << "integrate() done." << std::endl;
    std::cout << "steps_success = " << steps_success << std::endl;
    std::cout << "steps_fail = " << steps_fail << std::endl;
}

template<typename STEPPER>
bool ODD_integrator<STEPPER>::takeStep() {
    stepper->setBasePoint(t, x, f);
    bool step_success = false;
    // cut off stepsize to prevent overshoot
    if (t + stepper->next_h() > t_end) {
        stepper->next_h(t_end - t);
    }
    while (!step_success) {
        if (verbose > 0) std::cout << "\nnext stepsize: " << stepper->next_h() << std::endl;
        try {
            step_success = stepper->step();
        } catch (std::runtime_error e) {
            std::cerr << "ERROR: Caught Exception runtime_error:" << e.what() << std::endl;
            std::cerr << "stopping integration." << std::endl;
            break;
        }
        steps_fail += !step_success;
        steps_success += step_success;
        if (stepper->using_hmin()) {
            std::clog << "NOTE: Using minimal stepsize hmin = " << pars->find("hmin")->second << std::endl;
            break;
        }
    }
    if (!step_success) {
        std::cerr << "ERROR: point " << steps_success + 1 << " violates tolerance bounds.\nStopping integration..." << std::endl;
    } else {
        //        std::cout << "step successful. Accepting step..." << std::endl;
        // accept a successful step (add to solution)        
        acceptStep();
    }
    return step_success;

}


// integrate(t_end_new) can be used for extending the solution
// it continues the integration allow again maxsteps for the integration
// and resetting the step counters

template<typename STEPPER>
void ODD_integrator<STEPPER>::integrate(double t_end_new) {
    std::cout << "\n\nODD_integrator::integrate() called." << std::endl;
    if (t_end <= t) {
        std::cerr << "given t_end = " << t_end << " < " << t << " = actual t (nothing to do)." << std::endl;
        return;
    }
    t_end = t_end_new;
    steps_success = 0;
    steps_fail = 0;
    integrate();
}


// inspect_rhs() calls the user supplied right hand side and
// uses its footprint on xd to determine which delayed components
// are actually used. Only these will be stored and calculated.
//
// NOTE: To enable this check, I stripped the const-qualifier from
//       the xd parameter of the rhs function. This leaves room for the user
//       to change the given values 
//       (-> one should take care not to rely on the same interpolation results for two rhs calls)

template<typename STEPPER>
void ODD_integrator<STEPPER>::inspect_rhs() {
    if (verbose) std::cout << "inspect_rhs()" << std::endl;
    // clear xd to get the footprint
    xd->clear();
    // call rhs
    rhs(t, x, *xd, *rhs_pars, f);

    // trace which components were addressed
    if (verbose) {
        for (int j = 0; j != xd->size(); ++j) {
            for (auto p : (*xd)[j]) {
                std::cout << "rhs uses component xd[" << j << "][" << p.first << "]" << std::endl;
            }
        }
    }

}


// accept the last step taken by the stepper

template<typename STEPPER>
void ODD_integrator<STEPPER>::acceptStep() {
    t += stepper->last_h();
    x = stepper->lastPoint();

    // calculate new basepoint derivative
    sol->interpolate(t - *tau, *xd);
    xd->fillZeroDelays(x);
    rhs(t, x, *xd, *rhs_pars, f);

    // add point to soultion
    sol->push_back(t, x, f);
}


#endif	/* ODD_INTEGRATOR_H */

