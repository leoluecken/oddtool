/* 
 * File:   ODD_input.h
 * Author: luecken
 *
 * Created on 3. Februar 2014, 13:04
 */

#ifndef ODD_INPUT_H
#define	ODD_INPUT_H

#include <string>
#include <vector>
#include <utility>
#include <map>
#include <memory>
#include <boost/regex.hpp>
#include "ODD_utils.h"


/*
 * ODD_input 
 * 
 * parses configurations from the parameter files (default and user)
 * parses history data from the history file(s) if needed
 * 
 */

class ODD_input {
    
public:
    explicit ODD_input();    
    ~ODD_input();
    
    // init calls all load functions
    void init(std::string pfilename="ODD_parameters.txt", ODD_history_fct hist=nullptr);    
    
//    std::shared_ptr<const ODD_parameter> getParameters() const {return pars;}
    const ODD_parameter* getParameters() const {return pars.get();}
//    std::shared_ptr<const std::vector<double> > getRHSParameters() const {return rhs_pars;}
    const std::vector<double>* getRHSParameters() const {return rhs_pars.get();}
//    std::shared_ptr<const ODD_delays> getDelays() const {return tau;}
    const ODD_delays* getDelays() const {return tau.get();}
    const std::string& getName() const {return name;}
//    std::shared_ptr<const ODD_history> getHistory() const {return history;}
    const ODD_history* getHistory() const {return history.get();}
    const std::vector<size_t>* getUserDelayIndices() const {return user_delay_indices.get();}
    
    void printPars();

private:
    // load parameters from file parameter_filename
    void load_user_parameters();
    
    // load parameters from file default_parameter_filename
    void load_default_parameters();
    
    // load parameters in file fn into the map results
    void load_parameters(std::string fn, std::map<std::string,double>& results);
    
    // check if unknown parameters were specified and fill in missing pars with default values
    void check_and_complete_user_parameters();
    
    void check_and_init_rhs_pars(const std::map<size_t, double>& provided_pars);
    void check_and_init_tau(const std::map<size_t, double>& provided_lags);
    
    // load history from file history_filename
    void load_history();
    void load_history_from_single_file();
    void load_history_from_seperate_files();
        
    // name of the system
    std::string name="none";
    
    // track whether user specified a history for the derivatives
    bool derivatives_supplied=true;
    
    // parameters read from file parameter_filename
    std::unique_ptr<ODD_parameter> pars;
    
    // default parameters read from file default_parameter_filename 
    std::unique_ptr<ODD_parameter> default_parameters;
    
    // user parameters for the right hand side from file parameter_filename
    std::unique_ptr<std::vector<double> > rhs_pars;
    
    // list of delay values
    std::unique_ptr<ODD_delays> tau;
    
    // mapping that converts internal sorted indices to originally specified indices 
    std::unique_ptr<std::vector<size_t> > user_delay_indices;
    
    // history read from file 
    // first entry - timepoints
    // second entry - statepoints
    // third entry (if any) - derivative points     
    std::unique_ptr<ODD_history> history;    
    
    // file containing the default parameters data
    const std::string default_parameter_filename = "sources/ODD_default_parameters.txt";    
    
    // file containing the user parameters data
    std::string parameter_filename = "ODD_parameters.txt";
    
    // file(s) containing the initial data 
    std::string history_filename = "none";
    std::string history_t_filename = "none";
    std::string history_x_filename = "none";
    std::string history_f_filename = "none";
    
    
    // pattern to identify a numeric value, e.g. 100000, -3.5, +000023.23, 1.2e45... 
    const std::string numeric_value_pattern = "(?:-|\\+)?[[:d:]]+(?:\\.[[:d:]]*)?(e(?:-|\\+)[[:d:]]+)?";
    
    // pattern to identify a valid parameter definition
    const std::string c_comment = "/\\*[^/\\*]*\\*/";
    const std::string cpp_comment = "//[^/\\*]*";
    const std::string parameter_def_pattern = "[[:s:]]*([[:w:]]+)[[:s:]]*=[[:s:]]*([[:w:].]+|("+ numeric_value_pattern +"))[[:s:]]*;?[[:s:]]*(?:"+c_comment+"|"+cpp_comment+")?[[:s:]]*";
    
    // pattern to identify a delay name, e.g. t2 or tau_2 
    const std::string delay_name_pattern = "(?:t|tau)_?([[:d:]]+)";
    
    // pattern to identify an rhs parameter name, e.g. p1, p2, ... or p_1, p_2, ...
    const std::string rhs_parameter_name_pattern = "(?:p)[_]?([[:d:]]+)";
    
    // pattern to identify a good line in the history file 
    const std::string history_line_pattern = "([[:s:]]*"+numeric_value_pattern+"[[:s:]]*)+";
    
};


#endif	/* ODD_INPUT_H */

