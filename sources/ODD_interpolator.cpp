
#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include "ODD_interpolator.h"
#include "ODD_utils.h"

using namespace std;

typedef vector<double>::const_iterator c_iter;
typedef vector<double>::iterator iter;

ODD_interpolator::ODD_interpolator() : temp(4) {
}

void ODD_interpolator::init(size_t N) {
    this->N = N;
    hermite5_coeffs = vector<double>(N * 6);
}

/* This version of interpolate writes into x_query the requested interpolated values 
 * t_given, x_given are the known time and solution values, the derivatives in f_given
 * are used if the interpolation method_id is 'h' (Hermite interpolation)
 * t_query contains the values at which solution values are needed;
 * assuming t_given to be sorted strictly ascending
 * (TODO: use Newton interpolation)
 */
void ODD_interpolator::interpolate(
        c_iter t_given_begin, c_iter t_given_end, c_iter x_given_begin, c_iter f_given_begin,
        c_iter t_query_begin, c_iter t_query_end, ODD_delayed_values::iterator x_query_begin) {
//     cout << "called ODD_interpolator::interpolate() " << endl;

    // size of t_query
    unsigned long s = t_query_end - t_query_begin;

//    // trace
//    if (verbose>=0) {
//        cout << "size of query vector = " << s << endl;
//        cout << "t_query = ";
//        for (auto i = t_query_begin; i != t_query_end; ++i) {
//            cout << " " << *i;
//        }
//        cout << endl;
//        cout << "range of t_given = [" << *t_given_begin << "," << *(t_given_end - 1) << "]" << endl;
//    }

    double t; // temp var for actual timepoint to be interpolated
    auto j = t_given_begin; // temp var for iterator position

    for (size_t i = 0; i != s; ++i) {
        t = t_query_begin[i]; // actual time point for interpolation
//        if (verbose>=0) cout << "interpolating point " << t << endl;

        // use that t_query is strictly ascending
        if (j != t_given_begin) --j;
        j = find_interpolation_iterator(j, t_given_end, t);

        if (j == t_given_end) {
            // this return from find_interpolation_iterator() signals that t is not in the range of t_given => extrapolate
//            cerr << "NOTE: interpolation argument out of range (are you using h_max>tau_min?)." << endl;
            if(t<*t_given_begin) {
                j=t_given_begin;
            } else {
                j=t_given_end-1;
            }
            if (verbose>=1) cout << "extrapolating: " << t << " from " << *j << endl;
        } else {
            // point is in range => interpolate 
            if (verbose>=1) cout << "interpolating: next largest point in t_given: " << *j << endl;
        }
        
        interpolate_single(t, j, t_given_begin, t_given_end, x_given_begin, f_given_begin, x_query_begin + i);
    }
}

// interpolate (or extrapolate) a single value t with the given points around t_offset and nearest to t

void ODD_interpolator::interpolate_single(double t, c_iter t_offset, c_iter t_given_begin, c_iter t_given_end, c_iter x_given_begin, c_iter f_given_begin, ODD_delayed_values::iterator x_out) {

    bool extrapolate = (t < *t_given_begin || t >= *(t_given_end - 1));

    if (method_id == 'h') {
        if (extrapolate) {
            extrapolate_hermite5(t, t_offset, t_given_begin, t_given_end, x_given_begin, f_given_begin, x_out);
        } else {
            interpolate_hermite5(t, t_offset, t_given_begin, t_given_end, x_given_begin, f_given_begin, x_out);
        }
    } else if (method_id == 'l') {
        if (extrapolate) {
            extrapolate_linear(t, t_offset, t_given_begin, t_given_end, x_given_begin, x_out);
        } else {
            interpolate_linear(t, t_offset, t_given_begin, t_given_end, x_given_begin, x_out);
        }
    } else {
        throw invalid_argument("no valid id was supplied for the interpolation method.");
    }
}



// find iterator j corresponding to *(j-1) <= t < *j
// assuming a sorted list: *i < *j for i<j
// this method starts a recursive call of to its unchecked version,
// which doesn't test whether the given arguments are in range
// but conserves that *begin <= t < *end (see below)
// XXX: We could perhaps improve performance by estimating the bisection index

c_iter ODD_interpolator::find_interpolation_iterator(c_iter begin, c_iter end, double t) const {
    if (begin == end) return end;
    c_iter last = end - 1;
    if (t >= *last) {
        // failure
        // this includes the case begin==last (empty range)
        // if (verbose) cout << "t is out of range. Returning end." << endl;
        return end;
    } else if (t<*begin) {
        return begin;
    }
    
    if (last == begin + 1) return last; // base case -- success

    // bisect
    // note: we assured that (end - begin) >= 3
    auto it = begin + (end - begin) / 2;

    if (t >= *it) {
        return find_interpolation_iterator_unchecked(it, end, t);
    } else {
        return find_interpolation_iterator_unchecked(begin, it + 1, t);
    }
}


// unchecked version of bisection (never call with begin == end)
// assures for successive calls that t>=*begin and t<*end is preserved

c_iter ODD_interpolator::find_interpolation_iterator_unchecked(c_iter begin, c_iter end, double t) const {
    if (end == begin + 1) return end; // success
    auto it = begin + (end - begin) / 2; // bisect

    if (t >= *it) {
        return find_interpolation_iterator_unchecked(it, end, t);
    } else {
        return find_interpolation_iterator_unchecked(begin, it, t);
    }
}

void ODD_interpolator::interpolate_hermite5(double t, c_iter t_offset, c_iter t_given_begin, c_iter t_given_end, c_iter x_given_begin, c_iter f_given_begin, ODD_delayed_values::iterator x_out) {

//    if (verbose) cout << "ODD_interpolator::interpolate_hermite5() called" << endl;
    /* choose points for 3-point Hermite-interpolation (we know: t_offset>t_given_begin)
     * either (t_offset-2, t_offset-1, t_offset)
     * or (t_offset-1, t_offset, t_offset+1)                              
     * shift \in {0,1} selects the variant when using 
     * (t_offset-2, t_offset-1, t_offset) + shift */
    bool shift = 0;
    
    if (t_offset - 1 == t_given_begin) {
        // there's only one point to the left (which is t_offset-1)
        if (t_offset + 1 == t_given_end) {
            throw invalid_argument("Hermite interpolation needs at least 3 points.");
        }
        shift = 1;
    } else if (t_offset + 1 == t_given_end) {
        // there's only one point to the right (which is t_offset)
        shift = 0;
    } else {
        // determine shift such that |t - optional point| is minimized
        shift = (*(t_offset + 1) - t) < (t - *(t_offset - 2));
    }

//    if (verbose) cout << "shift = " << shift << endl;
//    cout << "*(t_offset-2) = " << *(t_offset-2) << endl;
//    cout << "*(t_offset-1) = " << *(t_offset-1) << endl;
//    cout << "*(t_offset) = " << *(t_offset) << endl;
//    cout << "*(t_offset+1) = " << *(t_offset+1) << endl; 

    double dt0 = *(t_offset + shift - 1)- *(t_offset + shift - 2);
    double dt1 = *(t_offset + shift)- *(t_offset + shift - 1);
    double dt2 = *(t_offset + shift)- *(t_offset + shift - 2);


    // next check is mainly to avoid division by zero,
    // ascending property is only checked on values, we work with
    if (dt0 <= 0 || dt1 <= 0 || dt2 <= 0) {
        throw invalid_argument("given time points are not ordered strictly ascending");
    }

    // make an index out of the iterator
    auto jx = t_offset - t_given_begin;

    // calculate coefficient of the interpolation polynom
    calculate_hermite5_coeffs(dt0, dt1, dt2, x_given_begin + (jx - 2 + shift) * N, f_given_begin + (jx - 2 + shift) * N, x_out);

    if (verbose>=2) print_hermite5_coeffs();

    // re-use dtj for storing offsets of grid points to interpolation point 
    dt0 = t - *(t_offset + shift - 2);
    dt1 = t - *(t_offset + shift - 1);
    dt2 = t - *(t_offset + shift);

    // evaluate interpolation polynomials and store results in x_out
    hermite5_poly(dt0, dt1, dt2, x_out);
}

void ODD_interpolator::extrapolate_hermite5(double t, c_iter t_offset, c_iter t_given_begin, c_iter t_given_end, c_iter x_given_begin, c_iter f_given_begin, ODD_delayed_values::iterator x_out) {

//    if (verbose) cout << "ODD_interpolator::extrapolate_hermite5() called" << endl;

    // determine shift {t_offset+shift-2, t_offset+shift-1, t_offset+shift} are the interpolation points
    int shift = 0;
    if (t_offset == t_given_begin) {
        shift = 2;
    } else if (t_offset == t_given_end - 1) {
        shift = 0;
    } else {
        throw range_error("extrapolate has to be called either with t_offset==t_given_begin or t_offset==*(t_given_end-1)");
    }

//    if (verbose) cout << "shift = " << shift << endl;

    double dt0 = *(t_offset + shift - 1)- *(t_offset + shift - 2);
    double dt1 = *(t_offset + shift)- *(t_offset + shift - 1);
    double dt2 = *(t_offset + shift)- *(t_offset + shift - 2);


    // next check is mainly to avoid division by zero,
    // ascending property is only checked on values, we work with
    if (dt0 <= 0 || dt1 <= 0 || dt2 <= 0) {
        throw invalid_argument("given time points are not ordered stirctly ascending");
    }

    // make an index out of the iterator
    auto jx = t_offset - t_given_begin;

    // calculate coefficient of the interpolation polynom
    calculate_hermite5_coeffs(dt0, dt1, dt2, x_given_begin + (jx - 2 + shift) * N, f_given_begin + (jx - 2 + shift) * N, x_out);

    if (verbose) print_hermite5_coeffs();

    // re-use dtj for storing offsets of grid points to interpolation point 
    dt0 = t - *(t_offset + shift - 2);
    dt1 = t - *(t_offset + shift - 1);
    dt2 = t - *(t_offset + shift);

    // evaluate interpolation polynomials and store results in x_out
    hermite5_poly(dt0, dt1, dt2, x_out);
}

/* evaluate the Hermite polynomials for the Ndim component functions. Given parameters:
 * dtj -- difference t-t_j of the interpolation point and the gridpoint t_j
 * x_out -- storage for the result
 * uses Hermite interpolation coefficients in hermite5_coeffs, which have to be 
 * calculated before by calling calculate_hermite5_coeffs()
 * (NOTE: the dtj arguments have a different meaning, there)
 */
void ODD_interpolator::hermite5_poly(double dt0, double dt1, double dt2, ODD_delayed_values::iterator x_out) {
//    if (verbose) cout << "hermite5_poly() called" << endl;
    auto& p = hermite5_coeffs;
    size_t i;
    for (auto& e : *x_out) {
        i = e.first;
        e.second = p[i * 6] + dt0 * (p[i * 6 + 1] + dt0 * (p[i * 6 + 2] + dt1 * (p[i * 6 + 3] + dt1 * (p[i * 6 + 4] + dt2 * p[i * 6 + 5]))));
    }
}

/* This calculates the coefficients of the fifth order approximation Hermite polynomials
 * for the Ndim component functions x_i at the points t0, t1, t2 provided the data  
 * dt0 = t1 - t0;  dt2 = t2 - t1;  dt2 = t2 - t0
 * x[j*Ndim + i] = x_i(tj)
 * f[j*Ndim + i] = x_i'(tj)
 * It uses generalized (Hermite) divided differences initializing the calculation at step P1.
 * The diagonal elements Pnn(i) are the coefficients of the i-th hermite polynomial.
 * Pnn(i) is stored in hermite5_coeffs[i * Ndim + n] for later use when evaluating the polynomial
 * x_out is a map whose keys indicate for which components the coefficients are needed
 */
void ODD_interpolator::calculate_hermite5_coeffs(double dt0, double dt1, double dt2, c_iter x, c_iter f, ODD_delayed_values::const_iterator x_out) {
//    if (verbose) cout << "calculate_hermite5_coeffs() called" << endl;
    //    for (int i = 0; i != Ndim; ++i) {
    size_t i;
//    if (verbose) cout << "number of components to interpolate: " << x_out->size() << endl;
    for (auto& e : *x_out) {
        i = e.first; // coordinate of the current component
//        if (verbose) cout << "calculating hermite coefficients for component " << i << "..." << endl;

        /* init coeffs of the i-th polynomial for state after first divided differences step P1 */
        hermite5_coeffs[i * 6] = x[i]; // == P00(i)
        hermite5_coeffs[i * 6 + 1] = f[i]; // == P11(i)
        hermite5_coeffs[i * 6 + 2] = (x[N + i] - x[i]) / dt0;
        hermite5_coeffs[i * 6 + 3] = f[N + i];
        hermite5_coeffs[i * 6 + 4] = (x[2 * N + i] - x[N + i]) / dt1;
        hermite5_coeffs[i * 6 + 5] = f[2 * N + i];

        /* calculate values for step P2 and store in temp */
        temp[0] = (hermite5_coeffs[i * 6 + 2] - hermite5_coeffs[i * 6 + 1]) / dt0;
        temp[1] = (hermite5_coeffs[i * 6 + 3] - hermite5_coeffs[i * 6 + 2]) / dt0;
        temp[2] = (hermite5_coeffs[i * 6 + 4] - hermite5_coeffs[i * 6 + 3]) / dt1;
        temp[3] = (hermite5_coeffs[i * 6 + 5] - hermite5_coeffs[i * 6 + 4]) / dt1;
        hermite5_coeffs[i * 6 + 2] = temp[0]; // == P22(i)

        /* calculate values for step P3 */
        hermite5_coeffs[i * 6 + 3] = (temp[1] - temp[0]) / dt0; // == P33(i)
        hermite5_coeffs[i * 6 + 4] = (temp[2] - temp[1]) / dt2;
        hermite5_coeffs[i * 6 + 5] = (temp[3] - temp[2]) / dt1;

        /* calculate values for step P4 and store in temp */
        temp[0] = (hermite5_coeffs[i * 6 + 4] - hermite5_coeffs[i * 6 + 3]) / dt2;
        temp[1] = (hermite5_coeffs[i * 6 + 5] - hermite5_coeffs[i * 6 + 4]) / dt2;
        hermite5_coeffs[i * 6 + 4] = temp[0]; // == P44(i)

        /* calculate last value for step P5 */
        hermite5_coeffs[i * 6 + 5] = (temp[1] - temp[0]) / dt2; // == P55(i)
    }
}

void ODD_interpolator::print_hermite5_coeffs() {
    print(hermite5_coeffs, 6);
}

/* write linearly interpolated values to x_out 
 * t_offset is the first position in ascending vector t_given such that *t_offset > t
 */
void ODD_interpolator::interpolate_linear(double t, c_iter t_offset, c_iter t_given_begin, c_iter t_given_end, c_iter x_given_begin, ODD_delayed_values::iterator x_out) {

    // next check is mainly to avoid division by zero
    // ascending property is only checked on values, we work with
    // XXX: we could drop this to improve performance
    if (*t_offset <= *(t_offset - 1)) {
        throw invalid_argument("given time points are not ordered stirctly ascending");
    }

    // make index of the iterator
    auto jx = t_offset - t_given_begin;

    // linear interpolation factors
    double q = (t - t_given_begin[jx - 1])*(t_given_begin[jx] - t_given_begin[jx - 1]);
    double p = 1 - q;

    // interpolate all components of x_query linearly (iterate through map<int, double> *x_out)
    size_t i;
    for (auto& e : *x_out) {
        i = e.first;
        e.second = p * x_given_begin[(jx - 1) * N + i] + q * x_given_begin[jx * N + i];
    }
}

/* write linearly extrapolated values to x_out */
void ODD_interpolator::extrapolate_linear(double t, c_iter t_offset, c_iter t_given_begin, c_iter t_given_end, c_iter x_given_begin, ODD_delayed_values::iterator x_out) {
    bool shift;
    if (t_offset == t_given_begin) {
        shift=1;
    } else if (t_offset == t_given_end - 1) {
        shift=0;
    } else {
        throw range_error("extrapolate has to be called either with t==*t_given_begin or t==*(t_given_end-1)");
    }

    // make index of the iterator
    auto jx = t_offset - t_given_begin + shift;

    // linear extrapolation factors
    double q = (t - t_given_begin[jx - 1])*(t_given_begin[jx] - t_given_begin[jx - 1]);
    double p = 1 - q;

    // extrapolate all components of x_query linearly (iterate through map<int, double> *x_out)
    size_t i;
    for (auto& e : *x_out) {
        i = e.first;
        e.second = p * x_given_begin[(jx - 1) * N + i] + q * x_given_begin[jx * N + i];
    }
}
