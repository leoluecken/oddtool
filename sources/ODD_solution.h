/* 
 * File:   ODD_solution.h
 * Author: luecken
 *
 * Created on 27. Januar 2014, 14:23
 * 
 * Solution object
 * holds a reference to the outfile, where the data is written to,
 * and a buffer containing the function history which might be needed to
 * interpolate delay values
 * Here and everywhere in the program, the buffer is constructed from 
 * two vectors for the time points and two vectors for the solution values
 * the solution value at time point t[i] is adressed by x[i*N + k], where 
 * k is the solution component and N is the system's dimension
 * 
 */

#ifndef ODD_SOLUTION_H
#define	ODD_SOLUTION_H

#include <vector>
#include <utility>
#include <map>
#include <string>
#include <memory>
#include "ODD_interpolator.h"
#include "ODD_output.h"
#include "ODD_utils.h"


typedef std::vector<std::vector<double> > ODD_history;
typedef std::map<std::string, double> ODD_parameter;

class SwapVectorStorage {
    double MINIMAL_RANGE; //MAX_DELAY; // time which has to be remembered for further evaluations of the rhs
    size_t MINIMAL_SAVE_POINTS = 100; // time which has to be remembered for further evaluations of the rhs
    size_t REQUIRED_INTERPOLATION_POINTS = 3;

    bool wi = 0; // write_index tells to which element of the storage pairs, 
    // (t_[0],t_[1]), (x_[0],x_[1]) and (f_[0],f_[1]) , new data is written, 
    // therefore write_index is either 0 or 1 (assured by type==bool)
    std::vector<std::vector<double> > t; // a pair of vectors of time points
    std::vector<std::vector<double> > x; // a pair of vectors (of N-dimensional points) for storing solution states
    std::vector<std::vector<double> > f; // a pair of vectors (of N-dimensional points) for storing values of the rhs
    std::vector<unsigned long> size; // a pair of write positions
    int N; // dimension of points in v_
    bool save_history;
    bool save_continuation_data;
    bool save_derivatives;
    double save_offset;
    struct {
        double offset;
        double period;
        double length;
        bool active;
        double last_offset;
    } partial_save;
    double save_step;
    double t_start;
    bool derivatives_available = false; // indicates, whether there are derivatives in the storage[!wi] which can be used for (hermite) interpolation
    bool using_history_fct = false;
    bool first_range = true; // indicates whether data[wi] is the first computed range (needed for discriminating interpolation case)
    bool last_range = false; // indicates the next flush is the last computed range before end (used in flush_all())
    bool saved_once = false; // indicates whether the first save is done
    double last_save_time; // invalid as long as saved_once==false

    std::string system_name;

    std::unique_ptr<ODD_output> out; // file output

    std::unique_ptr<ODD_interpolator> interpol; // interpolation

    ODD_history_fct* history_fct;

    /* connect() inits v_[write_index] with the last overlap points from data_[!write_index] 
     * and sets size_[write_index]=overlap */
    void connect(size_t overlap);
    void connect(); // uses overlap=REQUIRED_INTERPOLATION_POINTS    


    // flush contents of v_[!wi] if any and switch write_index
    void flush_and_swap();
    void flush(bool index);
    
    // flush a range in data[write_index] range
    void flush_range(c_iter t_begin, c_iter t_end, bool write_index, bool no_overlap = false);

    // save data for continuation of the solution to a single temp file
    void saveContinuationData();

    void initMembers(std::string name, const ODD_parameter* pars, size_t min_points);

    void initHistory(const ODD_history& history);

    void initHistory(double t0, std::vector<double> x0);

    void initHistory(ODD_history_fct history, double t_begin);
    
    void initPartialSave();

    // fills specified range of t and history(t) values into history
    void fillInHistory(double ts, double te, bool write_index, double h_step);

public:
    SwapVectorStorage();

    ~SwapVectorStorage();

    void init(std::string name, const ODD_parameter* pars, const ODD_history* history, size_t min_points = 10000);

    void init(std::string name, const ODD_parameter* pars, double t0, std::vector<double> x0, size_t min_points = 10000);

    void init(std::string name, const ODD_parameter* pars, ODD_history_fct hist, size_t min_points = 10000);

    unsigned long total_size() const {
        if (first_range) {
            return size[0] + size[1] - 1;
        } else {
            return size[0] + size[1] - REQUIRED_INTERPOLATION_POINTS;
        }
    }

    double operator[] (unsigned long j) const;

    void push_back(double t_new, const std::vector<double>& x_new, const std::vector<double>& f_new, bool swap_enabled = true);

    void interpolate(const std::vector<double>& t_query, ODD_delayed_values& x_query) const;

    // return the last x point
    std::vector<double> lastX() const;

    // return the last t point
    double lastT() const;

    size_t requiredInterpolationPoints() const {
        return REQUIRED_INTERPOLATION_POINTS;
    }

    // flush all contents of the storage
    void flush_all();

};

/* ODD_solution manages solution data contained in its data member */
class ODD_solution {
    size_t N; // the dimension of the system

    SwapVectorStorage data;

    void check_history(ODD_history hist);

public:
    ODD_solution();
    void init(std::string name, const ODD_parameter* pars, const ODD_history* history);
    void init(std::string name, const ODD_parameter* pars, ODD_history_fct history);

    // return the last x point 
    std::vector<double> lastX() const;

    // return the last t point
    double lastT() const;

    // get nr. of data points in data

    unsigned long size() const {
        return data.total_size();
    }

    size_t requiredInterpolationPoints() const {
        return data.requiredInterpolationPoints();
    }

    void push_back(double t_new, const std::vector<double>& x_new, const std::vector<double>& f_new);

    /* 
     * interpolate the solution history at a vector t_query of times 
     * writes interpolated solution values into x_query 
     */
    void interpolate(const std::vector<double>& t_query, ODD_delayed_values& x_query) const;
};

#endif	/* ODD_SOLUTION_H */

