/* 
 * File:   ODD_interpolator.h
 * Author: luecken
 *
 * Created on 7. Februar 2014, 11:41
 */



/* XXX:
 * (1) Interpolation can probably be optimized by using narrowing index ranges of all query times simultaneously 
 * and using monotonicity in some way to reduce comparisons
 * (2) provide version, which takes an initial guess for the ranges. This could speed up the repeated
 * query during the rk calculations
 * (3) using an estimate for the bisection point could reduce the bisection steps
 */


#ifndef ODD_INTERPOLATOR_H
#define	ODD_INTERPOLATOR_H

#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include "ODD_delayed_values.h"

typedef std::vector<double>::const_iterator c_iter;
typedef std::vector<double>::iterator iter;

class ODD_interpolator {
private:
    size_t N;
    std::vector<double> hermite5_coeffs; // storage for hermite coeffs of interpolation polynomials
    std::vector<double> temp; // temporary storage needed when calculating hermite5_coeffs

    // interpolation method 
    // h - hermite5
    // l - linear
    // TODO: n - newton (for init from discrete history data without specified derivatives)
    char method_id = 'l';

    // interpolate (or extrapolate) a single value t with the given points around t_offset and nearest to t
    void interpolate_single(double t, c_iter t_offset, c_iter t_given_begin, c_iter t_given_end, c_iter x_given_begin, c_iter f_given_begin, ODD_delayed_values::iterator x_out);

    // called only from find_interpolation_iterator() -- see below
    c_iter find_interpolation_iterator_unchecked(c_iter begin, c_iter end, double t) const;

    // write linearly interpolated values to x_out 
    // t_offset is the first position in ascending vector t_given such that *t_offset > t 
    void interpolate_linear(double t, c_iter t_offset, c_iter t_given_begin, c_iter t_given_end, c_iter x_given_begin, ODD_delayed_values::iterator x_out);
    // extrapolate linear
    void extrapolate_linear(double t, c_iter t_offset, c_iter t_given_begin, c_iter t_given_end, c_iter x_given_begin, ODD_delayed_values::iterator x_out);


    // Interpolate with degree 5 Hermite polynomial, using three values of x and three derivatives f(x).
    // Write results to x_out.
    // t_offset is the first position in ascending vector t_given such that *t_offset > t 
    void interpolate_hermite5(double t, c_iter t_offset, c_iter t_given_begin, c_iter t_given_end, c_iter x_given_begin, c_iter f_given_begin, ODD_delayed_values::iterator x_out);
    // Hermite_extrapolation
    void extrapolate_hermite5(double t, c_iter t_offset, c_iter t_given_begin, c_iter t_given_end, c_iter x_given_begin, c_iter f_given_begin, ODD_delayed_values::iterator x_out);

    // calculate coefficients for Hermite polynomial and store in hermite5_coeffs 
    void calculate_hermite5_coeffs(double dt0, double dt1, double dt2, c_iter x_begin, c_iter f_begin, ODD_delayed_values::const_iterator);

    // evaluate Hermite polynomials defined by hermite5_coeffs and write results to x_out
    void hermite5_poly(double dt0, double dt1, double dt2, ODD_delayed_values::iterator x_out);

    // print coeffs of hermite polynomial
    void print_hermite5_coeffs();

public:

    ODD_interpolator();
    void init(size_t N);
    size_t requiredInterpolationPoints() {return 3;}

    void setMethod(char id) {
        method_id = id;
    }

    // from the ascendingly sorted range [begin,end) return the first iterator i for a value v=*i > t
    c_iter find_interpolation_iterator(c_iter begin, c_iter end, double t) const;

    // interpolate values of x at t_query. (if f_given is provided, use information on derivatives)
    // results are written to x_query
    void interpolate(c_iter t_given_begin, c_iter t_given_end, c_iter x_given_begin, c_iter f_given_begin,
            c_iter t_query_begin, c_iter t_query_end, ODD_delayed_values::iterator x_query_begin);
};

#endif	/* ODD_INTERPOLATOR_H */

