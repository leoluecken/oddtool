/* 
 * File:   ODD_output.h
 * Author: luecken
 *
 * Created on 28. Januar 2014, 09:06
 * 
 * 
 */

#ifndef ODD_OUTPUT_H
#define	ODD_OUTPUT_H

#include <string>
#include <vector>
#include <fstream>
#include "ODD_utils.h"


/*
 * ODD_output
 * 
 * provides file-writing utilities, used by SwapVectorStorage
 * 
 */
class ODD_output {
    typedef std::vector<double>::const_iterator vec_iter;
public:
    ODD_output();
    void init(size_t N, bool save_derivatives=false, bool save_continuation_data=false, std::string filename_root="output");
    ~ODD_output();
        
    bool write(vec_iter t_begin, vec_iter t_end, vec_iter x_begin, vec_iter f_begin, int omit_last, int save_step = 0);
    bool write(vec_iter v_begin, vec_iter v_end, int data_dimension, std::string outfile_name, std::fstream::openmode m=std::fstream::out, bool continuation_data=true);
    bool write(vec_iter t_begin, vec_iter t_end, vec_iter x_begin, vec_iter f_begin, int data_dimension, std::string outfile_name, std::fstream::openmode m=std::fstream::out, bool continuation_data=true);
    bool write(vec_iter t_begin, vec_iter t_end, vec_iter x_begin, int data_dimension, std::string outfile_name, std::fstream::openmode m=std::fstream::out, bool continuation_data=true);
    bool write(vec_iter t_begin, vec_iter t_end, vec_iter x_begin, int omit_last=3, int save_step = 0);
    const std::string& getTempSubdirectory() const;
    const std::string& getFilenameRoot() const;
private:
    // filenaming
    std::string filename_root;
    std::string subdirectory;
    std::string continuation_data_subdirectory;
    std::string filename_suffix_x = "x";
    std::string filename_suffix_f = "f";
    std::string extension = ".txt";
    
    bool save_derivatives;
    
    // filestreams for timepoints, statepoints and derivativepoints
    std::ofstream fs_x;
    std::ofstream fs_f; 
    
    // system dimension
    int N;
};

#endif	/* ODD_OUTPUT_H */

