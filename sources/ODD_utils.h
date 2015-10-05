/* 
 * File:   ODD_utils.h
 * Author: luecken
 *
 * Created on 28. Januar 2014, 12:41
 * 
 * ODD_utils file defines some type aliases, vector arithmetic and some generic helper functions
 * 
 */


#ifndef ODD_UTILS_H
#define	ODD_UTILS_H

#include <vector>
#include <sstream>
#include <algorithm>
#include <limits>
#include <map>
#include <list>
#include <iostream>
#include <memory>
#include "ODD_delayed_values.h"


const static int verbose= 0; // control verbosity (so far verbose==0 means "not that much traces" and verbose!=0 means "all traces")
const double PI = 3.141592653589793238462;
const double EPS = std::numeric_limits<double>::epsilon();


class ODD_delayed_values;
typedef std::map<std::string, double> ODD_parameter;
typedef std::vector<std::vector<double> > ODD_history;
typedef std::vector<double> ODD_delays;
typedef void ODD_rhs(double t, const std::vector<double>& x, ODD_delayed_values& xd, const std::vector<double>& p, std::vector<double>& res);
typedef void ODD_history_fct(double t, std::vector<double>& res);

// print vector, n values in a row
void print(std::vector<double> v, std::size_t n);

// trace
void trace(std::string msg = "************** trace *************");



void wait();

// standard pattern to identify a line made of space-seperated numbers */
const std::string std_number_line_pattern = "([[:s:]]*(?:-|\\+)?[[:d:]]+(?:\\.[[:d:]]*)?(e(?:-|\\+)[[:d:]]+)?[[:s:]]*)+";

// load all numbers in file filename, provided a pattern to recognize a line of numbers
std::vector<double> load_all_numbers(std::string filename, std::string number_line_pattern = std_number_line_pattern);


// vector output

template<typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> v) {
    os << "[";
    for (auto& t : v) os << t << " ";
    os << "]";
    return os;
}

// vector arithmetic
std::vector<double> operator+(double x, const std::vector<double>& v);
std::vector<double> operator+(const std::vector<double>& v, double x);
std::vector<double> operator-(double x, const std::vector<double>& v);
std::vector<double> operator-(const std::vector<double>& v, double x);
std::vector<double> operator*(double x, const std::vector<double>& v);
std::vector<double> operator*(const std::vector<double>& v, double x);
std::vector<double> operator+(const std::vector<double>& v, const std::vector<double>& w);
std::vector<double> operator-(const std::vector<double>& v, const std::vector<double>& w);

// componentwise absolut values
std::vector<double> cwise_abs(const std::vector<double>& v);
// componentwise division
std::vector<double> cwise_divide(const std::vector<double>& v, const std::vector<double>& w);


// get a vector of the keynames from a map

template<typename KEY, typename VAL>
std::vector<KEY> getKeys(const std::map<KEY, VAL>& m) {
    std::vector<KEY> res(m.size());
    for (auto p : m) {
        res.push_back(p.first);
    }
    return res;
}


// check whether an element is in a container

template<typename CONTAINER, typename E>
bool is_in(const CONTAINER& c, const E& e) {
    return std::find(c.begin(), c.end(), e) != c.end();
}

// parse one word from the source value to the target value via stringstream

template<typename SOURCE_TYPE, typename TARGET_TYPE>
TARGET_TYPE& parse(const SOURCE_TYPE& s, TARGET_TYPE& t) {
    std::stringstream ss;
    ss << s;
    ss >> t;
    return t;
}

template<typename A, typename B>
std::pair<B, A> flip_pair(const std::pair<A, B> &p) {
    return std::pair<B, A > (p.second, p.first);
}

template<typename A, typename B>
std::map<B, A> flip_map(const std::map<A, B> &src) {
    std::map<B, A> dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()),
            flip_pair<A, B>);
    return dst;
}

/* determines number of different elements in v*/
template<typename T>
size_t nr_of_distinct_elements(const std::vector<T> v) {
    std::map<T, size_t> m;
    for (T e : v) {
        m[e] += 1;
    }
    return m.size();
}

/* determines number of different elements in m*/
template<typename A, typename B>
size_t nr_of_distinct_elements(const std::map<A, B> m) {
    std::map<B, size_t> n;
    for (auto e : m) {
        n[e.second] += 1;
    }
    return n.size();
}

#endif	/* ODD_UTILS_H */

