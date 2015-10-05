/* 
 * File:   ODD_stepper.h
 * Author: luecken
 *
 * Created on 28. Januar 2014, 09:05
 * 
 * The stepper mainly provides the function step(),
 * which carries out an integration step, handles step sizes
 * and provides error estimates
 * 
 */

#ifndef ODD_STEPPER_H
#define	ODD_STEPPER_H

#include <vector>
#include <memory>
#include "ODD_utils.h"
#include "ODD_solution.h"

class ODD_stepper {
public:
    ODD_stepper();
    virtual ~ODD_stepper() {}
    virtual void init(ODD_rhs user_rhs, const ODD_parameter* pars, const ODD_solution* sol, const ODD_delays* tau, const std::vector<double>* p)=0;
    
    virtual bool step() = 0;
    virtual bool step(double h) = 0;
    virtual double next_h() const {return h_next;}
    virtual double last_h() const {return h_last;}
    virtual void next_h(double h) { h_next=h;}
    virtual void setHmax(double h) { hmax=h;}
    virtual double getHmax() const { return hmax;}
    virtual const std::vector<double>& lastPoint() const {return last_point;}
    bool using_hmin() const {return h_last == hmin;}
    
protected:
    size_t N; // system dimension
    
    // relative and absolute tolerances
    double rel_tol, abs_tol; 
    
    // stepwidth of the last calculation
    double h_last;
    // proposed stepwidth for the next calculation
    double h_next;
    // minimal and maximal stepwidth
    double hmin, hmax;

    /* value of the last calculated step*/
    std::vector<double> last_point;

    ODD_rhs* user_rhs;
    
    /* user supplied rhs */
    virtual void rhs(double t, const std::vector<double>& x, ODD_delayed_values& xd, const std::vector<double>& p, std::vector<double>& result) {
        // switch to user indices while calling user supplied rhs
        xd.useUserIndices(true);
        user_rhs(t, x, xd, p, result);
        xd.useUserIndices(false);
    }
    
    // pointer to solution owned by the integrator (used to obtain delayed values by interpolation)
    const ODD_solution* sol;
    
    // pointer to the delay times owned by the integrator
    const ODD_delays* tau; 
    // indicates from which index i on the delays are 0
    size_t tau0_offset; 
    
    // pointer to the user_specified rhs-parameters owned by the integrator
    const std::vector<double>* p; 
};



struct CKButcherTableau {    
    // ----- method parameters (Cash-Karp) -----
    const double a2 = 1.0 / 5.0;
    const double a3 = 3.0 / 10.0;
    const double a4 = 3.0 / 5.0;
    const double a5 = 1.0;
    const double a6 = 7.0 / 8.0;

    const double b21 = 1.0 / 5.0;
    const double b31 = 3.0 / 40.0;
    const double b32 = 9.0 / 40.0;
    const double b41 = 3.0 / 10.0;
    const double b42 = -9.0 / 10.0;
    const double b43 = 6.0 / 5.0;
    const double b51 = -11.0 / 54.0;
    const double b52 = 5.0 / 2.0;
    const double b53 = -70.0 / 27.0;
    const double b54 = 35.0 / 27.0;
    const double b61 = 1631.0 / 55296.0;
    const double b62 = 175.0 / 512.0;
    const double b63 = 575.0 / 13824.0;
    const double b64 = 44275.0 / 110592.0;
    const double b65 = 253.0 / 4096.0;

    const double c41 = 37.0 / 378.0;
    const double c43 = 250.0 / 621.0;
    const double c44 = 125.0 / 594.0;
    const double c46 = 512.0 / 1771.0;
    const double c51 = 2825.0 / 27648.0;
    const double c53 = 18575.0 / 48384.0;
    const double c54 = 13525.0 / 55296.0;
    const double c55 = 277.0 / 14336.0;
    const double c56 = 1.0 / 4.0;
};



class CashKarpStepper : public ODD_stepper {
public:
    explicit CashKarpStepper();
    ~CashKarpStepper() {}
    void init(ODD_rhs user_rhs, const ODD_parameter* pars, const ODD_solution* sol, const ODD_delays* tau, const std::vector<double>* p) override;
    void init_xd(ODD_delayed_values* xd) {this->xd=xd;};
    void setBasePoint(double t, const std::vector<double>& x, const std::vector<double>& f);
    bool step() override;
    bool step(double new_h) override;
//    void rhs(double t, const std::vector<double>& x, const std::vector<double>& xd, const std::vector<double>& p, std::vector<double>& result) override;
private:
    double safety_factor=0.9;
    double minimal_stepsize_rescale=0.1;
    double maximal_stepsize_rescale=10.0;
    
    /* base time t0 */
    double t0;
    /* base point x0 */
    std::vector<double> x0;
    /* base value f(x0) */
    std::vector<double> f0;
    /* CashKarp Butcher Tableau */
    CKButcherTableau bt;
    /* temp values for runge kutta calculations */
    std::vector<double> k1, k2, k3, k4, k5, k6, x_tmp;
    
    /* temp container for retrieval of delayed solution values */
    ODD_delayed_values* xd;
    
    /* temp containers for increments */
    std::vector<double> increment_order_4, increment_order_5;
    /* temp_containers for componentwise errors and scales */
    std::vector<double> abs_err, rel_err, component_scales;
};
typedef CashKarpStepper CKStepper;




#endif	/* ODD_STEPPER_H */

