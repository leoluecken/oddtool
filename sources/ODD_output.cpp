#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
//#include <cstdlib> // system

#include "ODD_output.h"
#include "ODD_utils.h"

using namespace std;

typedef std::vector<double>::const_iterator vec_iter;

ODD_output::ODD_output() {
}


// init() makes subdirectories and opens filestreams for solution data flushes
// filestreams for continuation data are opened and closed each time they are saved 
// (which should be fewer times)

void ODD_output::init(size_t N, bool save_deriv, bool save_continuation_data, string name) {
    // TODO: check whether files exists. In case print warning and ask for override permission 

    cout << "ODD_output::init()\nfilename root = " << filename_root << endl;

    this->N = N;
    save_derivatives = save_deriv;
    filename_root = name;
    subdirectory = name + "_output/"; // holds the complete output
    continuation_data_subdirectory = subdirectory + "continuation_data/"; // holds continuation data

    // create subdirectories
    system(("mkdir " + subdirectory).c_str());

    if (save_continuation_data) {
        system(("mkdir " + continuation_data_subdirectory).c_str());
    }

    // set full precision for writing doubles
    fs_x.precision(std::numeric_limits<double>::digits10);
    if (save_derivatives) {
        fs_f.precision(std::numeric_limits<double>::digits10);
    }
    // open files
    cout << "opening filestreams ... " << endl;
    //    cout << "dimension: " << N << endl;
    fs_x.open(subdirectory + filename_suffix_x + extension, ofstream::out);
    if (save_derivatives) {
        fs_f.open(subdirectory + filename_suffix_f + extension, ofstream::out);
    }
}

// destructor closes filestreams

ODD_output::~ODD_output() {
    // close filestreams
    fs_x.close();
    if (save_derivatives) {
        fs_f.close();
    }
    //    cout << "~ODD_output()" << endl;
}

const string& ODD_output::getTempSubdirectory() const {
    return continuation_data_subdirectory;
}

const string& ODD_output::getFilenameRoot() const {
    return filename_root;
}

// write() writes data from range [v_begin,v_end) in data_dimension columns to the outfile

bool ODD_output::write(vec_iter v_begin, vec_iter v_end, int data_dimension, std::string outfile_name, fstream::openmode m, bool continuation_data) {
    if (continuation_data) outfile_name = continuation_data_subdirectory + outfile_name;
    //    cout << "ODD_output::write() called (to file " << outfile_name << ")" << endl;
    ofstream fs;
    fs.open(outfile_name, m);
    fs.precision(std::numeric_limits<double>::digits10);
    int count = 0;
    while (v_begin != v_end) {
        fs << fixed << *(v_begin++) << " ";
        ++count;
        //        cout<<count<<endl;
        if (count % data_dimension == 0) fs << endl;
    }
    fs.close();
    return fs.good();
}

// writes data from ranges [t_begin,t_end), [x_begin,x_begin+t_end-t_begin), and [f_begin,f_begin+t_end-t_begin)
// in 1 + 2*data_dimension columns to the outfile
// requires data_dimension > 0

bool ODD_output::write(vec_iter t_begin, vec_iter t_end, vec_iter x_begin, vec_iter f_begin, int data_dimension, std::string outfile_name, std::fstream::openmode m, bool continuation_data) {
    if (continuation_data) outfile_name = continuation_data_subdirectory + outfile_name;
    //    cout << "ODD_output::write() called (to file " << outfile_name << ")" << endl;
    ofstream fs;
    fs.open(outfile_name, m);
    fs.precision(std::numeric_limits<double>::digits10);
    int count = 0;
    while (t_begin != t_end) {
        fs << fixed << *(t_begin++) << " ";
        for (count = 0; count != data_dimension; ++count) {
            fs << fixed << *(x_begin++) << " ";
        }
        for (count = 0; count != data_dimension; ++count) {
            fs << fixed << *(f_begin++) << " ";
        }
        fs << endl;
    }
    fs.close();
    return fs.good();
}

// writes data from ranges [t_begin,t_end) and [x_begin,x_begin+t_end-t_begin)
// in 1 + data_dimension columns to the outfile requires data_dimension > 0

bool ODD_output::write(vec_iter t_begin, vec_iter t_end, vec_iter x_begin, int data_dimension, std::string outfile_name, std::fstream::openmode m, bool continuation_data) {
    if (continuation_data) outfile_name = continuation_data_subdirectory + outfile_name;
    //    cout << "ODD_output::write() called (to file " << outfile_name << ")" << endl;
    ofstream fs;
    fs.open(outfile_name, m);
    fs.precision(std::numeric_limits<double>::digits10);
    int count = 0;
    while (t_begin != t_end) {
        fs << fixed << *(t_begin++) << " ";
        for (count = 0; count != data_dimension; ++count) {
            fs << fixed << *(x_begin++) << " ";
        }
        fs << endl;
    }
    fs.close();
    return fs.good();
}


// writes time points in range [t_begin, t_end) and corresponding x's and f's to 
// the corresponding filestreams; Is called repeatedly (at every swap event) from SwapVectorStorage

bool ODD_output::write(vec_iter t_begin, vec_iter t_end, vec_iter x_begin, vec_iter f_begin, int omit, int save_step /* = 0 */) {
    //    cout << "ODD_output::write() called (to files " << subdirectory << "*.txt)" << endl;
    if (save_derivatives) {
        // write derivative data values
        int fcount = 0;
        auto t_iter = t_begin;
        double t_last = *t_iter - save_step; // sth <= tbegin-save_step to trigger saving of the first point
        //    for (unsigned long j = 0; j != size - write_last; ++j) {
        while (t_iter != t_end - omit) {
            if (t_last <= *t_iter - save_step) {
                t_last = *t_iter;
                // save
                fs_f << fixed << *t_iter++ << " ";
                for (fcount = 0; fcount != N; ++fcount) {
                    fs_f << fixed << *f_begin++ << " ";
                }
                fs_f << endl;
            }else {
            ++t_iter;
            f_begin+=N;
        }
        }
    }
    // write x and t data
    return write(t_begin, t_end, x_begin, omit, save_step);
}


// writes time points in range [t_begin, t_end) and corresponding x's and f's to 
// the corresponding filestreams; Is called repeatedly (at every swap event) from SwapVectorStorage

bool ODD_output::write(vec_iter t_begin, vec_iter t_end, vec_iter x_begin, int omit, int save_step /* = 0 */) {
//    cout << "ODD_output::write() called (to files " << subdirectory << "*.txt) with save_step = "<< save_step << endl;


    // write data values
    int xcount = 0;
    double t_last = *t_begin - save_step; // sth <= tbegin-save_step to trigger saving of the first point
    //    for (unsigned long j = 0; j != size - write_last; ++j) {
    while (t_begin != t_end - omit) {
        if (t_last <= *t_begin - save_step) {
            t_last = *t_begin;
            fs_x << fixed << *t_begin++ << " ";
            for (xcount = 0; xcount != N; ++xcount) {
                //            cout << *x_begin << " ";
                fs_x << fixed << *x_begin++ << " ";
            }
            fs_x << endl;
            //        cout << *t_begin << endl;
        } else {
            ++t_begin;
            x_begin+=N;
        }
    }
    return fs_x.good();
}

