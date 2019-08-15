
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <utility>
#include <memory>
#include "ODD_solution.h"
#include "ODD_output.h"
#include "ODD_utils.h"

using namespace std;

SwapVectorStorage::SwapVectorStorage() : t(2), x(2), f(2), size(2) {
}

void SwapVectorStorage::init(string name, const ODD_parameter* pars, const ODD_history* history, size_t min_points) {
    cout << "SwapVectorStorage::init()" << endl;
    //    if (verbose) cout << "creating SwapVectorStorage " << filename << "..." << endl;
    //    if (verbose) cout << "dimension: " << pars->find("N")->second << endl;

    initMembers(name, pars, min_points);
    initHistory(*history);
    initPartialSave();
}

void SwapVectorStorage::init(string name, const ODD_parameter* pars, double t0, vector<double> x0, size_t min_points) {
    cout << "SwapVectorStorage::init()" << endl;
    //    if (verbose) cout << "creating SwapVectorStorage " << filename << "..." << endl;
    //    if (verbose) cout << "dimension: " << pars->find("N")->second << endl;

    initMembers(name, pars, min_points);
    initHistory(t0, x0);
    initPartialSave();
}

void SwapVectorStorage::init(string name, const ODD_parameter* pars, ODD_history_fct history, size_t min_points) {
    cout << "SwapVectorStorage::init()" << endl;
    //    if (verbose) cout << "creating SwapVectorStorage " << filename << "..." << endl;
    //    if (verbose) cout << "dimension: " << pars->find("N")->second << endl;

    initMembers(name, pars, min_points);
    initHistory(history, pars->find("t_start")->second);
    if (verbose >= 0)cout << "Using user supplied history function." << endl;
    initPartialSave();
}

// init partial save configuration (call immediately after initHistory() to assure t0= t[wi][0])

void SwapVectorStorage::initPartialSave() {
    double t0 = max(t[wi][0], save_offset); // start time 
    double det = partial_save.offset - t0;
    partial_save.active = partial_save.period > 0;
    partial_save.last_offset = t0 + fmod(det, partial_save.period) + (det > 0 ? -partial_save.period : 0);
}

void SwapVectorStorage::initMembers(string name, const ODD_parameter* pars, size_t min_points) {
    N = pars->find("N")->second;
    system_name = name;
    save_derivatives = pars->find("save_derivatives")->second;
    save_continuation_data = pars->find("save_continuation_data")->second;
    save_history = pars->find("save_history")->second;
    save_step = pars->find("save_step")->second;
    save_offset = pars->find("save_offset")->second;
    t_start = pars->find("t_start")->second;

    partial_save.offset = pars->find("partial_save_offset")->second;
    partial_save.period = pars->find("partial_save_period")->second;
    partial_save.length = pars->find("partial_save_length")->second;

    //    interpol = make_shared<ODD_interpolator > ();
    interpol = unique_ptr<ODD_interpolator > (new ODD_interpolator());
    interpol->init(N);
    REQUIRED_INTERPOLATION_POINTS = interpol->requiredInterpolationPoints();

    //    out = make_shared<ODD_output > ();
    out = unique_ptr<ODD_output > (new ODD_output());
    out->init(N, pars->find("save_derivatives")->second, pars->find("save_continuation_data")->second, name);

    // set parameters for storage management
    MINIMAL_RANGE = pars->find("tau_max")->second + pars->find("hmax")->second * 2;
    MINIMAL_SAVE_POINTS = min_points;
    N = pars->find("N")->second;

    //    if (verbose) cout << "reserving storage..." << endl;

    // reserve minimal length for the storage vectors and init first elements (presupposed by connect())
    t[0].reserve(MINIMAL_SAVE_POINTS);
    t[0] = {0};
    t[1].reserve(MINIMAL_SAVE_POINTS);
    t[1] = {0};
    x[0].reserve(N * MINIMAL_SAVE_POINTS);
    x[0] = vector<double>(N);
    x[1].reserve(N * MINIMAL_SAVE_POINTS);
    x[1] = vector<double>(N);
    f[0].reserve(N * MINIMAL_SAVE_POINTS);
    f[0] = vector<double>(N);
    f[1].reserve(N * MINIMAL_SAVE_POINTS);
    f[1] = vector<double>(N);
    //    cout << "reserved.\nsetting write position ..." << endl;
    // set write positions
    size[0] = 0;
    size[1] = 0;
    // direct write operations to v_[0] initially
    wi = 0;
}

// TODO: make sure |data|>=3

void SwapVectorStorage::initHistory(const ODD_history& history) {
    t[!wi] = history[0];
    x[!wi] = history[1];
    if (history.size() == 3) {
        // history of derivatives was supplied
        f[!wi] = history[2];
        derivatives_available = true;
    } else {
        derivatives_available = false;
    }
    size[!wi] = history[0].size();
    connect(1);
    cout << "initial size of storage: " << this->total_size() << endl;
}


// init with constant initial data

void SwapVectorStorage::initHistory(double t0, vector<double> x0) {
    vector<double> f0(N);
    for (int n = 0; n != REQUIRED_INTERPOLATION_POINTS; ++n) {
        push_back(t0 - MINIMAL_RANGE + MINIMAL_RANGE * n / ((double) REQUIRED_INTERPOLATION_POINTS - 1), x0, f0, false);
    }
    derivatives_available = true;
    wi = !wi;
    connect(1);
    cout << "initial size of storage: " << this->total_size() << endl;
}

// init with user supplied history function

void SwapVectorStorage::initHistory(ODD_history_fct history, double t_begin) {
    history_fct = history;
    using_history_fct = true;
    derivatives_available = false;
    if (save_history) {
        // put some values of the history in the storage, if they are of interest.
        // they wont be used in calculations, 
        fillInHistory(t_begin - MINIMAL_RANGE, t_begin, wi, min(save_step, MINIMAL_RANGE / (double) REQUIRED_INTERPOLATION_POINTS));
    } else {
        // only push initial point
        fillInHistory(t_begin, t_begin, wi, 1);
    }
    wi = !wi;
    connect(1);
    cout << "initial size of storage: " << this->total_size() << endl;
}

// fills specified range of t and history(t) values into history

void SwapVectorStorage::fillInHistory(double ts, double te, bool write_index, double h_step) {
    //    if(verbose>=0) cout<<"fillInHistory("<<ts<<", "<<te<<")"<<endl; 
    if (h_step <= 0) throw invalid_argument("fillInHistory: h_step has to be positive.");
    if (ts == te) clog << "NOTE: called fillInHistory with one point range." << endl;
    bool wi_orig = wi;
    wi = write_index; // direct pushes to write_index
    vector<double> xs(N);
    vector<double> xf(N);
    while (abs(te - ts) >= EPS) {
        history_fct(ts, xs);
        push_back(ts, xs, xf, false);
        ts = min(te, ts + h_step);
    }
    history_fct(te, xs);
    push_back(te, xs, xf, false);
    wi = wi_orig; // reset wi to original position
}

SwapVectorStorage::~SwapVectorStorage() {
    //    if (verbose >= 0) cout << "entering ~SwapVectorStorage()" << endl;

    //    if (verbose >= 0) cout << "storage sizes:\nsize[0] = " << size[0] << "\nsize[1] = " << size[1] << endl;
    //    if (verbose >= 0) cout << "storage ranges:\nrange[0] = [" << t[0][0] << ", " << t[0][size[0] - 1] << "]\nrange[1] = [" << t[1][0] << ", " << t[1][size[1] - 1] << "]" << endl;

    // save continuation data
    if (save_continuation_data) saveContinuationData();

    // write all data remaining in the storage to output file
    flush_all();
    //    if (verbose >= 0) cout << "leaving ~SwapVectorStorage()" << endl;
}


// save data in the last time interval of length MINIMAL_RANGE

void SwapVectorStorage::saveContinuationData() {
    //    if (verbose) cout << "save_continuation_data()" << endl;
    bool save_derivatives = true; // this overwrites save_derivatives specified by user for his output

    // minimal time point to be saved
    auto tmin = t[wi][size[wi] - 1] - MINIMAL_RANGE;
    //    if (verbose) cout << "required tmin = " << tmin << endl;
    //    if (verbose) cout << "given tmin = " << t[!wi][0] << endl;

    if (tmin < t_start) {
        clog << "WARNING: no continuation data will be saved since integration time was less than MIN_RANGE." << endl;
        return;
    }
    // check if data is split over both storage vectors
    bool split = (tmin < t[wi][0]);
    vector<double>::const_iterator it = t[wi].begin();
    fstream::openmode write_mode = fstream::out;
    size_t j;
    if (split) {
        if (using_history_fct) {
            //            if (verbose) cout << "collecting continuation data from provided history function with grid_size = " << save_step << endl;
            // XXX: no information on derivatives of the provided history fct is available, 
            //I make no attempt to approximate, since this case here is not very important 
            save_derivatives = false;
            clog << "WARNING not writing derivative continuation data. Integration was shorter than SwapVectorStorage::MINIMAL_RANGE." << endl;
            return;
        } else {
            // if (verbose) cout << "collecting continuation data from data[!write_index]" << endl;
            // find tmin in t_[!write_index]
            it = interpol->find_interpolation_iterator(t[!wi].begin(), t[!wi].begin() + size[!wi], tmin);
            // iterator corresponding to the last point before that value
            if (it != t[!wi].begin()) --it;
        }
        // corresponding index j
        j = it - t[!wi].begin();
        // write data points in data_[!write_index] to files
        size_t overlap = first_range ? 1 : REQUIRED_INTERPOLATION_POINTS;
        if (save_derivatives) {
            out->write(t[!wi].begin() + j, t[!wi].begin() + size[!wi] - overlap, x[!wi].begin() + N*j, f[!wi].begin() + N*j, N, system_name + "_history.txt", write_mode);
        } else {
            out->write(t[!wi].begin() + j, t[!wi].begin() + size[!wi] - overlap, x[!wi].begin() + N*j, N, system_name + "_history.txt", write_mode);
        }
        write_mode = fstream::app;
    }
    // set it to the right position in t_[write_index]
    if (!split) {
        if (verbose) cout << "collecting continuation data only from data[write_index]" << endl;
        if (verbose) cout << "looking for tmin = " << tmin << " in [" << *(t[wi].begin()) << ", " << *(t[wi].begin() + size[wi] - 1) << "]" << endl;
        // find tmin in t_[write_index]
        it = interpol->find_interpolation_iterator(t[wi].begin(), t[wi].begin() + size[wi], tmin);
        if (verbose) cout << "found next largest entry at t = " << *it << ", that is, at index j = " << it - t[wi].begin() << endl;
        // iterator corresponding to the last point before that value
        if (it != t[wi].begin()) --it;

    } else {
        if (verbose) cout << "collecting continuation data from data[write_index]" << endl;
        // don't save the 1st value (it is contained int data[!write_index] and was saved)
        it = t[wi].begin() + 1;
    }
    // corresponding index j
    j = it - t[wi].begin();
    // write data points in data[write_index] to files
    if (save_derivatives) {
        out->write(t[wi].begin() + j, t[wi].begin() + size[wi], x[wi].begin() + N*j, f[wi].begin() + N*j, N, system_name + "_history.txt", write_mode);
    } else {
        out->write(t[wi].begin() + j, t[wi].begin() + size[wi], x[wi].begin() + N*j, N, system_name + "_history.txt", write_mode);
    }
}


// [j] returns j-th time value

double SwapVectorStorage::operator[] (unsigned long j) const {
    if (j >= size[0] + size[1] - 1 || j < 0) {
        throw invalid_argument("SwapVectorStorage: index out of bounds");
    }
    if (j < size[!wi]) {
        return t[!wi][j];
    } else {
        return t[wi][j - size[!wi]];
    }
}

void SwapVectorStorage::flush_and_swap() {
    if (verbose >= 1) cout << "entering flush_and_swap() at t = " << t[wi][size[wi] - 1] << endl;
    //        if (verbose >= 0) cout << "size[write_index=" << wi << "] = " << size[wi] << endl;
    //    if (verbose >= 2) cout << "size[!write_index=" << !wi << "] = " << size[!wi] << endl;
    //    if (verbose >= 2) cout << "x[" << !wi << "].capacity() = " << x[!wi].capacity() << endl;
    // flush    
    flush(!wi);
    // swap
    wi = !wi;
    connect();
}


// flush the specified index

void SwapVectorStorage::flush(bool index) {
    if (size[index] == 0) {
        return;
    }
    if (t[index][size[index] - 1] > save_offset) {
        c_iter t_begin = t[index].begin();
        if (!saved_once) {
            // first time saving: find offset
            t_begin = interpol->find_interpolation_iterator(t_begin, t_begin + size[index], save_offset);
            if (t_begin != t[index].begin()) --t_begin;
        }

        if (partial_save.active) {
            cout << "flushing partial data..." << endl;
            double next_off = partial_save.last_offset + partial_save.length;
            if (next_off < *t_begin) next_off += partial_save.period;
            double next_on = partial_save.last_offset + partial_save.period;


            auto t_end = t[index].begin() + size[index];
            bool no_overlap = false;
            while (t_begin != t_end) {
                cout << "t = " << *t_begin << endl;
                cout << "next off-event at t = " << next_off << endl;
                cout << "next on-event at t = " << next_on << endl;
                if (next_off < next_on) {
                    // saving is on... find position of the next off-event
                    c_iter j = interpol->find_interpolation_iterator(t_begin, t_end, next_off);
                    if (j != t_end) {
                        ++j;
                        no_overlap = true;
                    }
                    // ... and flush
                    cout << "saving range [" << *t_begin << ", " << *(j - 1) << "]" << endl;
                    flush_range(t_begin, j, index, no_overlap);
                    t_begin = j;
                    no_overlap = false;
                    if (t_begin != t_end) next_off += partial_save.period;
                } else {
                    // saving is off. just advance t_begin until min(t_end, next_on)
                    t_begin = interpol->find_interpolation_iterator(t_begin, t_end, next_on);
                    if (t_begin != t_end) {
                        --t_begin;
                        partial_save.last_offset = next_on;
                        next_on += partial_save.period;
                    }
                }
            }
        } else {
            // no partial saving - just flush the whole storage
            auto t_end = t[index].begin() + size[index];
            flush_range(t_begin, t_end, index);
        }

        saved_once = true;
    } else {
//        cout << "skip output: t[index][end] < save_offset" << endl;
    }
    size[index] = 0;
}

// flush data[write_index] in range [*t_begin, *(t_end-1)]

void SwapVectorStorage::flush_range(c_iter t_begin, c_iter t_end, bool write_index, bool no_overlap /* = false */) {
    size_t overlap = (no_overlap || last_range) ? 0 : (first_range ? 1 : REQUIRED_INTERPOLATION_POINTS);

    // index of first save point
    auto i_begin = t_begin - t[write_index].begin();

    // flush data[write_index]
    // save the new write vector and reset its write position
    bool success;
    if (save_derivatives) {
        if (!derivatives_available) {
            char msg[200];
            sprintf(msg, "WARNING: derivatives can not be saved because not available in [%f,%f]", t[write_index][0], t[write_index][size[write_index] - 1]);
            throw runtime_error(msg);
        }
        success = out->write(t_begin, t_end, x[write_index].begin() + N*i_begin, f[write_index].begin() + N*i_begin, overlap, save_step);
    } else {
        success = out->write(t_begin, t_end, x[write_index].begin() + N*i_begin, overlap, save_step);
    }
    if (!success) throw runtime_error("could not evoke ODD_output::write successfully");
}

void SwapVectorStorage::flush_all() {
    cout << "entering flush_all() at t = " << t[wi][size[wi] - 1] << endl;
    // flush trailing storage
    flush(!wi);
    // flush leading storage
    flush(wi);
}

/* reset data[write_index] by copying the last overlap elements of data[!write_index] 
 * to first elements of data[write_index] */
void SwapVectorStorage::connect() {
    connect(REQUIRED_INTERPOLATION_POINTS);
}

void SwapVectorStorage::connect(size_t overlap) {
    if (verbose > 0) cout << "connecting..." << endl;
    //    if (verbose) cout << "copying last elements..." << endl;
    // copy last history elements
    size_t s = size[!wi];
    if (overlap > s) throw invalid_argument("connect() requires overlap <= size[!wi]");
    for (size_t n = overlap; n != 0; --n) {
        push_back(t[!wi][s - n],{x[!wi].begin()+(s - n) * N, x[!wi].begin()+(s - n + 1) * N},
        {
            f[!wi].begin()+(s - n) * N, f[!wi].begin()+(s - n + 1) * N
        });
    }
    //    if (verbose >= 0) cout << "t_[wi][0 : 2] = [" << t[wi][0] << ", " << t[wi][1] << ", " << t[wi][2] << "]" << endl;
}

void SwapVectorStorage::push_back(double t_new, const std::vector<double>& x_new, const std::vector<double>& f_new, bool swap_enabled) {
    if (verbose > 0) cout << "adding " << size[wi] + 1 << "-th point to SwapVectorStorage::data[" << wi << "]..." << endl;
    //    if (verbose == 0) cout << "t = " << t_new << "\nx = [" << x_new[0] << ",..., " << x_new[N - 1] << "]" << "\nf = [" << f_new[0] << ",..., " << f_new[N - 1] << "]" << endl;
    if (verbose > 0) cout << "t = " << t_new << "\nx = [" << x_new[0] << ", " << x_new[1] << ", " << x_new[2] << "]" << "\nf = [" << f_new[0] << ", " << f_new[1] << ", " << f_new[2] << "]" << endl;
    // add point
    if (t[wi].size() > size[wi]) {
        t[wi][size[wi]] = t_new;
        for (size_t k = 0; k != N; ++k) {
            x[wi][N * size[wi] + k] = x_new[k];
            f[wi][N * size[wi] + k] = f_new[k];
        }
    } else {
        t[wi].push_back(t_new);
        for (size_t k = 0; k != N; ++k) {
            x[wi].push_back(x_new[k]);
            f[wi].push_back(f_new[k]);
            //            if(verbose) cout << "pushed x = " << x[wi][x[wi].size()-1] <<endl;
        }
    }
    ++size[wi];
    //    if (verbose) cout << "   new size = " << total_size() << "\n   size_[write_index] = " << size[wi] << endl;
    //        if (verbose) cout << "t_[write_index][last]-t_[write_index][0] = " << t[wi][size[wi] - 1] - t[wi][0] << endl;

    if (!swap_enabled) {
        return;
    }

    // reorganize storage if appropriate
    if ((size[wi] >= MINIMAL_SAVE_POINTS) && t_new - t[wi][0] > MINIMAL_RANGE) {
        flush_and_swap();
        derivatives_available = true;
        first_range = false;
        using_history_fct = false; // from first save on, the user supplied hist_fct is not used anymore
    }

}

vector<double> SwapVectorStorage::lastX() const {
    //    if(verbose>0) cout<<"lastX() of data["<<wi<<"] with size[wi] = "<<size[wi]<<endl;
    //    cout << "x[wi][0] = " <<x[wi][0]<<endl;
    //    cout << "x[!wi][0] = " <<x[!wi][0]<<endl;
    return
    {
        x[wi].begin() + N * (size[wi] - 1), x[wi].begin() + N * size[wi]
    };
}

double SwapVectorStorage::lastT() const {
    return t[wi][size[wi] - 1];
}

ODD_solution::ODD_solution() {
}

void ODD_solution::init(std::string name, const ODD_parameter* pars, const ODD_history* history) {
    //    cout << "ODD_solution::init()" << endl;
    N = pars->find("N")->second;
    size_t min_save_points = (size_t) pars->find("min_points")->second;
    check_history(*history);
    if ((*history)[0].size() == 1) {
        // const initial data
        data.init(name, pars, (*history)[0][0], (*history)[1], min_save_points);
    } else {
        data.init(name, pars, history, min_save_points);
    }
}

void ODD_solution::init(std::string name, const ODD_parameter* pars, ODD_history_fct history) {
    //    cout << "ODD_solution::init()" << endl;
    N = pars->find("N")->second;
    size_t min_save_points = (size_t) pars->find("min_points")->second;
    data.init(name, pars, history, min_save_points);
}

void ODD_solution::check_history(ODD_history hist) {
    bool derivatives_supplied = (hist.size() == 3);
    //    if (verbose) {
    //        cout << "history[0].size() = " << hist[0].size() << endl;
    //        cout << "history[1].size() = " << hist[1].size() << endl;
    //        if (derivatives_supplied) {
    //            cout << "history[2].size() = " << hist[2].size() << endl;
    //        }
    //    }

    if (hist[0].size() * N != hist[1].size()) {
        stringstream ss("error in history data: N*<nr of timepoints> (=");
        ss << N * hist[0].size() << ") is not equal to <nr of state points>/N (=" << hist[1].size() << ")";
        throw invalid_argument(ss.str());
    }

    if (derivatives_supplied && (hist[0].size() * N != hist[2].size())) {
        stringstream ss("error in history data: N*<nr of timepoints> (=");
        ss << N * hist[0].size() << ") is not equal to <nr of derivative points>/N (=" << hist[1].size() << ")";
        throw invalid_argument(ss.str());
    }
}

vector<double> ODD_solution::lastX() const {
    return data.lastX();
}

double ODD_solution::lastT() const {
    return data.lastT();
}

void ODD_solution::push_back(double t_new, const std::vector<double>& x_new, const std::vector<double>& f_new) {
    // add point
    data.push_back(t_new, x_new, f_new);
}

// this interpolates the values x_query based on t_ and x_
// assumes strictly ascending query times t_query (this corresponds to a corresponding
// descending ordering of the delays which should be assured at system initialization)
// TODO: provide another argument which specifies the components for which the interpolated 
// values are needed (or change x_query type and include it there)

void ODD_solution::interpolate(const vector<double>& t_query, ODD_delayed_values& x_query) const {
    //    if (verbose) {
    //        cout << "ODD_solution::interpolate() called with t_query =";
    //        for (auto t : t_query) {
    //            cout << " " << t;
    //        }
    //        cout << endl;
    //    }
    if (t_query.size() != 0) data.interpolate(t_query, x_query);
}



// this interpolates the values x_query based on t_ and x_
// assumes strictly ascending query times t_query (this corresponds to a corresponding
// descending ordering of the delays which should be assured at system initialization)

void SwapVectorStorage::interpolate(const vector<double>& t_query, ODD_delayed_values& x_query) const {
    //    cout << t_query.size() << endl;
    if (verbose > 0) cout << "SwapVectorStorage::interpolate() with ";
    //        cout  << "t_[write_index][0] = " << t[wi][0] << endl;
    if (verbose > 0) cout << "t_query[0] = " << t_query[0] << ", t_query[1] = " << t_query[1] << endl;

    // find position s that splices the query into query from data_[write_index] and data_[!write_index] 
    auto s = t_query.begin();
    auto t_query_end = t_query.begin() + x_query.getTau0_offset(); // ignore zero-delays during interpolation

    if (t[wi][0] > *t_query.begin()) {
        s = interpol->find_interpolation_iterator(t_query.begin(), t_query_end, t[wi][0]);
    }
    // if the point *(s-1) is actually the connecting point, it should be searched in the second part
    if (s != t_query.begin() && *(s - 1) == t[wi][0]) --s;

    //     if(s!=t_query.begin()) {
    //         cout << "seperating iterator refers to one after " << *(s - 1) << endl;
    //     } else {
    //         cout << "seperating iterator refers to " << *s << endl;
    //     }


    if (s == t_query_end) {
        //        if (verbose) cout << "interpolation case: all in data_[!write_index]" << endl;
        // all is before t_[write_index]
        if (using_history_fct) {
            // no interpolation, just call the history fct
            vector<double> v(N); // tmp storage
            // iterator through x_query.data
            auto j = x_query.begin();
            for (auto i = t_query.begin(); i != s; ++i) {
                // get history values at *i = t-tau_i
                history_fct(*i, v);
                // visit all needed delayed components x[i][k]=x_k(t-t_i) for time *i = t-t_i 
                // and store the calculated value (ck.first==index, ck.second==value)
                for (auto& ck : *j) {
                    ck.second = v[ck.first];
                }
                // step forward in x_query.data
                ++j;
            }
        } else {
            if (derivatives_available) {
                interpol->setMethod('h');
            } else {
                interpol->setMethod('l');
            }
            // interpolation using the stored timepoints
            interpol->interpolate(t[!wi].begin(), t[!wi].begin() + size[!wi],
                    x[!wi].begin(), f[!wi].begin(),
                    t_query.begin(), t_query_end, x_query.begin());
        }
    } else if (s == t_query.begin()) {
        //        if (verbose) cout << "interpolation case: all in data_[write_index]" << endl;
        // all is in t_[write_index] (or larger)
        interpol->setMethod('h');
        interpol->interpolate(t[wi].begin(), t[wi].begin() + size[wi],
                x[wi].begin(), f[wi].begin(),
                t_query.begin(), t_query_end, x_query.begin());
    } else {
        //        if (verbose) cout << "interpolation case: general, t_query distributed in t_[0][:] and t_[1][:]" << endl;
        // general situation: two calls to interpolate()        
        // first interpolate the part which lies in data[!write_index] 
        // return value j points to the first position of the second part
        //        cout << "interpolating part < " << *s << " in range [" << *(t_[!write_index].begin()) << ", " << *(t_[!write_index].begin() + size_[!write_index] - 1) << ")" << endl;
        if (using_history_fct) {
            // no interpolation, just call the history fct
            vector<double> v(N); // tmp storage
            // iterator through x_query.data
            auto j = x_query.begin();
            for (auto i = t_query.begin(); i != s; ++i) {
                // get history values at *i = t-tau_i
                history_fct(*i, v);
                // visit all needed delayed components x[i][k]=x_k(t-t_i) for time *i = t-t_i 
                // and store the calculated value (ck.first==index, ck.second==value)
                for (auto& ck : *j) {
                    ck.second = v[ck.first];
                }
                // step forward in x_query.data
                ++j;
            }
        } else {
            if (derivatives_available) {
                interpol->setMethod('h');
            } else {
                interpol->setMethod('l');
            }
            interpol->interpolate(t[!wi].begin(), t[!wi].begin() + size[!wi],
                    x[!wi].begin(), f[!wi].begin(), t_query.begin(), s, x_query.begin());
        }
        // then interpolate the second part (which lies in data[write_index])
        //        cout << "interpolated first part." << endl;
        //        cout << "interpolating part >= " << *s << " in range [" << *(t_[write_index].begin()) << ", " << *(t_[write_index].begin() + size_[write_index] - 1) << ")" << endl;
        interpol->setMethod('h');
        interpol->interpolate(t[wi].begin(), t[wi].begin() + size[wi],
                x[wi].begin(), f[wi].begin(), s, t_query_end, x_query.begin() + (s - t_query.begin()));
    }

    if (verbose > 0) x_query.print();

}

