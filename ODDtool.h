/* 
 * File:   ODDtool.h
 * Author: luecken
 *
 * Created on 18. März 2014, 14:02
 */

#ifndef ODDTOOL_H
#define	ODDTOOL_H

#include "sources/ODD_integrator.h"

/* signature of the user supplied specification of the right hand side of \dot{x}= f(t,x,xd;p)
 * t  - actual time
 * x  - actual state of variables
 * xd - delayed states of variables (xd[i][j] is the value of the j-th component at time t-tau_i)
 * p  - user specified parameters 
 * NOTE: the value of f(t,x,xd;p) should be written to the N-dimensional vector result
 * NOTE: xd should not supposed to be changed though it is not read-only. 
 *       The program exploits this to determine dynamically which delayed components are actually needed.
 */
void rhs_function(double t, const std::vector<double>& x, ODD_delayed_values& xd, const std::vector<double>& p, std::vector<double>& result);

/* signature of the user-supplied history function 
 * result is an N-dimensional vector to which history_function(t, result) should write
 * the state of the system at time t <= t_start
 */
void history_function(double t, std::vector<double>& result);


#endif	/* ODDTOOL_H */

