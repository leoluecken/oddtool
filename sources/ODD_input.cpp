#include <string>
#include <sstream>
#include <exception>
#include <vector>
#include <utility>
#include <map>
#include <boost/regex.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <memory>
#include "ODD_delayed_values.h"
#include "ODD_input.h"

using namespace std;
using namespace boost;

ODD_input::ODD_input() {
}

void ODD_input::init(string pfilename, ODD_history_fct hist) {
    parameter_filename = pfilename;

    // init parameter and history pointers
//    pars = make_shared<ODD_parameter > ();
    pars = unique_ptr<ODD_parameter>(new ODD_parameter());
//    default_parameters = make_shared<ODD_parameter > ();
    default_parameters = unique_ptr<ODD_parameter>(new ODD_parameter());
//    rhs_pars = make_shared<vector<double> >();
    rhs_pars = unique_ptr<vector<double> >(new vector<double>());
//    tau = make_shared<ODD_delays > ();
    tau = unique_ptr<ODD_delays>(new ODD_delays());
//    user_delay_indices = make_shared<vector<size_t> >();
    user_delay_indices = unique_ptr<vector<size_t> >(new vector<size_t>());
//    history = make_shared<ODD_history > ();
    history = unique_ptr<ODD_history>(new ODD_history());

    // load parameters
    load_default_parameters();
    load_user_parameters();
    check_and_complete_user_parameters();

    if (verbose) printPars();

    if (hist == nullptr || pars->find("force_history_from_file")->second != 0) {
        // load history from file(s)
        load_history();
    } else {
        cout << "user supplied history function." << endl;
    }
}

// dtor: no work since working with smart ptrs

ODD_input::~ODD_input() {
//    cout << "~ODD_input()" << endl;
}

// load parameters from file parameter_filename 

void ODD_input::load_parameters(string fn, ODD_parameter& results) {
    clog << "loading parameters from: " << fn << endl;

    // open file fn
    ifstream fs;
    fs.open(fn);
    if (!fs.good()) throw runtime_error("Could not open parameter file: " + fn);

    string line;
    double d_buffer;
    int i_buffer;

    // user specified rhs pars 
    map<size_t, double> provided_rhs_pars;

    // user specified delays 
    map<size_t, double> provided_lags;

    // init regex's
    regex r(parameter_def_pattern), r2(delay_name_pattern), r3(rhs_parameter_name_pattern), r4("history(_[[:alpha:]])?_filename");
    smatch match, match2;

    // parse parameters in the file
    while (getline(fs, line)) {
        if (regex_match(line, match, r)) {
            /* found a parameter definition */
            string parameter_name = match.str(1);
            string parameter_value = match.str(2);
            if (regex_match(parameter_name, match2, r4)) {
                /* user specified history file(s) */
                string history_suffix = match2.str(1);
                if (history_suffix == "") {
                    history_filename = parameter_value;
                } else if (history_suffix == "_t") {
                    history_t_filename = parameter_value;
                } else if (history_suffix == "_x") {
                    history_x_filename = parameter_value;
                } else if (history_suffix == "_f") {
                    history_f_filename = parameter_value;
                }
            } else if (parameter_name == "name") {
                name = parameter_value;
            } else if (regex_match(parameter_name, match2, r2)) {
                /* this is a delay specification */
                int tau_nr = parse(match2.str(1), i_buffer);
                double tau_val = parse(parameter_value, d_buffer);
                provided_lags.insert({tau_nr, tau_val});
            } else if (regex_match(parameter_name, match2, r3)) {
                /* this is a user parameter specification */
                int j = parse(match2.str(1), i_buffer);
                provided_rhs_pars[j] = parse(parameter_value, d_buffer);
            } else {
                /* all other allowed parameters have double values */
                results[parameter_name] = parse(parameter_value, d_buffer);
                cout << "loaded pars[\"" << parameter_name << "\"] = " << results[parameter_name] << endl;
            }
        }
    }

    if (!is_in(getKeys(results), "N")) {
        throw runtime_error("please provide system dimension N in parameter file.");
    }

    if (provided_lags.size() != 0) {
        check_and_init_tau(provided_lags);
    } else {
        clog << "no delay specifications in file " << fn << endl;
        results["tau0_offset"]= 0;
    }

    // check provided_rhs_pars and copy to member rhs_pars
    if (provided_rhs_pars.size() != 0) {
        check_and_init_rhs_pars(provided_rhs_pars);
    } else {
        clog << "no rhs parameter specifications in file " << fn << endl;
    }
}


// init delays: sort descendingly and merge equal delays into one index

void ODD_input::check_and_init_tau(const map<size_t, double>& provided_lags) {
    clog << "check_and_init_tau()" << endl;
    size_t Nt = provided_lags.size();
    if (Nt == 0) {
        clog << "no delays provided." << endl;
        return;
    }
    auto& tau_ref = *tau;
    auto& pars_ref = *pars;

    // merge equal user provided delays into descendingly sorted list tau and remember user indices 
    // (to access delayed values in function call with the user provided function)
    sort_and_make_index_map(provided_lags, tau_ref, *user_delay_indices);
    
    // determine first delay index at which delay is not positive (deal with delays tau_i==0)
    size_t j=0;
    while(tau_ref[j]>0 && j!= tau_ref.size()) ++j;
    pars_ref["tau0_offset"]= j;

    clog << "nr of delays provided: " << Nt << endl;
    clog << "nr of distinct delays: " << tau_ref.size() << endl;
    if (tau_ref.size() != 0 && *(tau_ref.rbegin()) < 0) {
        stringstream ss;
        ss << "All delays have to be non-negative. (tau_min = " << *(tau_ref.rbegin()) << ")";
        throw range_error(ss.str());
    }
    clog << "nr of zero delays: " << tau_ref.size()-j << endl;

    // get minimal and maximal positive delay
    pars_ref["tau_max"] = *(tau_ref.begin());
    auto tmin_it = --tau_ref.end();
    while (*tmin_it == 0 && tmin_it != tau_ref.begin()) --tmin_it;
    pars_ref["tau_min"] = *tmin_it;
    clog << "tau_min = " << pars_ref["tau_min"] << endl;
    clog << "tau_max = " << pars_ref["tau_max"] << endl;
}


// put provided rhs parameters in a vector; if some index within maximal and minimal 
// specified indices remains unspecified, it is set to 0.

void ODD_input::check_and_init_rhs_pars(const map<size_t, double>& provided_rhs_pars) {
    clog << "check_and_init_rhs_pars()" << endl;
    clog << "nr of rhs parameters provided: " << provided_rhs_pars.size() << endl;
    if (provided_rhs_pars.size() == 0) {
        return;
    }
    // collect provided indices 
    vector<size_t> provided_rhs_pars_keys = getKeys(provided_rhs_pars);
    // nr of parameters == largest given index + 1 
    int np = *max_element(provided_rhs_pars_keys.begin(), provided_rhs_pars_keys.end()) + 1;
    // clog << "maximal index provided: " << (np-1) << endl;
//    rhs_pars = make_shared<vector<double> >(np);
    rhs_pars = unique_ptr<vector<double> >(new vector<double>(np));
    for (int j = 0; j != np; ++j) {
        if (is_in(provided_rhs_pars_keys, j)) {
            (*rhs_pars)[j]= provided_rhs_pars.find(j)->second;
            cout << "p[" << j << "] = " << provided_rhs_pars.find(j)->second << endl;
        } else {
            (*rhs_pars)[j] = 0;
            clog << "WARNING: rhs parameter nr " << j << " was not specified.\nDefault initializing to 0" << endl;
        }
    }
}

void ODD_input::load_default_parameters() {
    load_parameters(default_parameter_filename, *default_parameters);
}

void ODD_input::load_user_parameters() {
    load_parameters(parameter_filename, *pars);
}


// check if unknown parameters were specified and fill in missing pars with default values

void ODD_input::check_and_complete_user_parameters() {
    cout << "check_and_complete_user_parameters()" << endl;
    // list of parameters (obtained from default_parameters) 
    vector<string> par_keys = getKeys(*default_parameters);
    
    // list of user-specified parameters
    vector<string> given_keys = getKeys(*pars);

    // check which parameters were provided by the user
    // warn which specified parameters are unknown 
    for (auto k : given_keys) {
        if (!is_in(par_keys, k)) {
            clog << "WARNING: user specification of unknown parameter: " << k << endl;
        }
    }

    // default init for pars which are uninitialized
    for (auto k : par_keys) {
        if (!is_in(given_keys, k)) {
            clog << "parameter " << k << " was not specified by user\ninitialization to default value " << (*default_parameters)[k] << endl;
            (*pars)[k] = (*default_parameters)[k];
        }
    }

    /* check condition for interpolating delays from history 
     * TODO: handle violations by allowing extrapolation in interpolation algorithm
     */
    if ((*pars)["hmax"] >= (*pars)["tau_min"]) {
        clog << "WARNING: hmax > tau_min..." << endl;
    }
    
    // check save_step
    if((*pars)["save_step"] < 0){
        throw invalid_argument("save_step < 0");
    }
    
    // check consistency of save_offset and save_history
    if ((*pars)["save_history"] == 0 && (*pars)["save_offset"] < (*pars)["t_start"]){
        clog << "WARNING: supplied save_offset < t_start, but save_history == 0. Setting save_offset = t_start."<<endl;
        (*pars)["save_offset"] = (*pars)["t_start"];
    }
    
    // check partial save config
    if((*pars)["partial_save_period"] < 0){
        throw invalid_argument("partial_save_period < 0");
    }    
    if((*pars)["partial_save_period"] > 0 && (*pars)["partial_save_length"] == 0 ){
        clog << "WARNING: automatically disabling partial save because partial_save_length == 0" <<endl;
        (*pars)["partial_save_period"]=0;
    }
    if((*pars)["partial_save_period"] > 0 && (*pars)["partial_save_length"] > (*pars)["partial_save_period"]){
        clog << "WARNING: automatically disabling partial save because partial_save_length > partial_save_period" <<endl;
        (*pars)["partial_save_period"]=0;
    }
    if((*pars)["partial_save_period"] > 0 && (*pars)["partial_save_length"] == (*pars)["partial_save_period"]){
        clog << "WARNING: automatically disabling partial save because partial_save_length == partial_save_period" <<endl;
        (*pars)["partial_save_period"]=0;
    }
}

void ODD_input::load_history() {
    // check whether filenames were supplied by user
    bool load_from_multi_files = history_t_filename != "none" && history_x_filename != "none";

    if (load_from_multi_files) {
        clog << "loading history from:\n" << history_t_filename << "\n" << history_x_filename << "\n" << history_f_filename << endl;
        load_history_from_seperate_files();
    } else {
        clog << "loading history from: " << history_filename << endl;
        load_history_from_single_file();
    }
}

// load history from files history_<id>_filename 
// and write data into history

void ODD_input::load_history_from_seperate_files() {
    derivatives_supplied = (history_f_filename != "");

    vector<double> t0 = load_all_numbers(history_t_filename, history_line_pattern);
    vector<double> x0 = load_all_numbers(history_x_filename, history_line_pattern);
    vector<double> f0;
    if (derivatives_supplied) {
        f0 = load_all_numbers(history_f_filename, history_line_pattern);
    }

//    cout << "size t0 = " << t0.size() << endl;

    history->push_back(t0);
    history->push_back(x0);
    if (derivatives_supplied) {
        history->push_back(f0);
    }

    bool sizes_match = (t0.size()*(*pars)["N"] == x0.size());
    if (derivatives_supplied) {
        sizes_match = sizes_match && (f0.size() == x0.size());
    }

    if (!sizes_match) {
        /* no matches in the whole file */
        clog << "WARNING: inconsistent data lengths:\n" <<
                history_t_filename << ": " << t0.size() << " points.\n" <<
                history_x_filename << ": " << x0.size() << " points." << endl;
        if (derivatives_supplied) {
            clog << history_f_filename << ": " << f0.size() << " points." << endl;
        }
    } else {
        clog << "Found " << t0.size() << " history data points." << endl;
        clog << "derivatives_supplied = " << derivatives_supplied << endl;
        if (t0[t0.size() - 1] - t0[0] < pars->find("tau_max")->second) {
            clog << "WARNING: provided history time is smaller than maximal delay tau_max = "<< pars->find("tau_max")->second <<"." << endl;
        }
    }
}

// load history from file history_filename 
// and write data into history

void ODD_input::load_history_from_single_file() {

    ifstream fs;
    fs.open(history_filename);
    if (!fs.good()) {
        clog << "Could not open history file: " << history_filename << endl;
        return;
    }

    vector<double> t0;
    vector<double> x0;
    vector<double> f0;

    string line;
    stringstream ss;
    double buffer;
    long good_line_count = 0;

    regex r(history_line_pattern);
    bool match;
    size_t N = pars->find("N")->second;

    while (getline(fs, line)) {
        // check whether this line is a valid history line
        try {
            match = regex_match(line, r);
        } catch (std::exception& e) {
            clog << "catching exception thrown by regex_match():\n" << e.what() << endl;
            match = false;
        }
        if (!match) {
            continue;
        }
        /* process line */
        ss.clear();
        ss << line;
        ss >> buffer;
        t0.push_back(buffer);
        for (int j = 0; j < N; ++j) {
            if (ss >> buffer) {
                x0.push_back(buffer);
            } else {
                throw invalid_argument("history read-in error: state point dimension is not equal to N = " + N);
            }
        }
        ++good_line_count;
        if (derivatives_supplied) {
            // d._h._s. is initialized = true
            // we enter at least once to check whether derivatives are supplied
            for (int j = 0; j < N; ++j) {
                if (ss >> buffer) {
                    f0.push_back(buffer);
                } else {
                    derivatives_supplied = false;
                    break;
                }
            }
        }
    }
    fs.close();

    history->push_back(t0);
    history->push_back(x0);
    if (derivatives_supplied) {
        history->push_back(f0);
    }

    if (good_line_count == 0) {
        /* no matches in the whole file */
        clog << "WARNING: Could not find any history data in file: " + history_filename << endl;
    } else {
        clog << "Found " << good_line_count << " lines of history data." << endl;
        clog << "derivatives_supplied = " << derivatives_supplied << endl;
        if (t0[t0.size() - 1] - t0[0] < pars->find("tau_max")->second) {
            clog << "WARNING: provided history time (" << t0[0] << "," << t0[t0.size() - 1] << ") is smaller than maximal delay tau_max = "<< pars->find("tau_max")->second <<"." << endl;
        }
    }
}

void ODD_input::printPars() {
    cout << "\nparameters given:" << endl;
    for (auto p : *pars) {
        cout << "pars[\"" << p.first << "\"] = " << p.second << "\n";
    }
    cout << endl;
}