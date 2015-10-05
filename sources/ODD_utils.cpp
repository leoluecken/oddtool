
#include <vector>
#include <utility>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <boost/regex.hpp>
#include <iostream>
#include <cmath>

using namespace std;
using namespace boost;

// print vector, n values in a row

void print(std::vector<double> v, std::size_t n) {
    for (size_t i = 0; i < v.size(); ++i) {
        if (i % n == 0) std::cout << "\n";
        std::cout << std::fixed << v[i] << " ";
    }
    std::cout << "\n" << endl;
}

void trace(string msg = "************** trace *************") {
    std::cout << msg << std::endl;
}



void wait(){
    string s;
    cout<<"press enter..."<<endl;
    getline(cin,s);
}

// load all numbers in file filename, provided a pattern to recognize a line of numbers

vector<double> load_all_numbers(string filename, string number_line_pattern) {
    ifstream fs;
    string line;
    stringstream ss;
    double buffer;
    regex r(number_line_pattern);
    bool match;
    vector<double> res;

    fs.open(filename);
    if (!fs.good()) {
        throw runtime_error("Could not open file: " + filename);
    }
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
        while (ss >> buffer) {
            // cout<<"found value "<<buffer<< endl;
            res.push_back(buffer);
        }
    }
    fs.close();
    return res;
}


// vector arithmetic

vector<double> operator+(double x, const vector<double>& v) {
    int n = v.size();
    vector<double> ret;
    ret.reserve(n);
    for (unsigned long i = 0; i != n; ++i) {
        ret.push_back(v[i] + x);
    }
    return ret;
}

vector<double> operator+(const vector<double>& v, double x) {
    return x + v;
}

vector<double> operator-(double x, const vector<double>& v) {
    int n = v.size();
    vector<double> ret;
    ret.reserve(n);
    for (unsigned long i = 0; i != n; ++i) {
        ret.push_back(x - v[i]);
    }
    return ret;
}

vector<double> operator-(const vector<double>& v, double x) {
    return -x + v;
}

vector<double> operator*(double x, const vector<double>& v) {
    int n = v.size();
    vector<double> ret;
    ret.reserve(n);
    for (unsigned long i = 0; i != n; ++i) {
        ret.push_back(v[i] * x);
    }
    return ret;
}

vector<double> operator*(const vector<double>& v, double x) {
    return x*v;
}

vector<double> operator+(const vector<double>& v, const vector<double>& w) {
    int n = v.size();
    int m = w.size();
    vector<double> ret;
    ret.reserve(n);
    for (unsigned long i = 0; i != n; ++i) {
        ret.push_back(v[i] + w[i]);
    }
    return ret;
}

vector<double> operator-(const vector<double>& v, const vector<double>& w) {
    int n = v.size();
    vector<double> ret;
    ret.reserve(n);
    for (unsigned long i = 0; i != n; ++i) {
        ret.push_back(v[i] - w[i]);
    }
    return ret;
}


// componentwise absolut values

vector<double> cwise_abs(const vector<double>& v) {
    int n = v.size();
    vector<double> ret;
    ret.reserve(n);
    for (unsigned long i = 0; i != n; ++i) {
        ret.push_back(fabs(v[i]));
    }
    return ret;
}


// componentwise division

vector<double> cwise_divide(const vector<double>& v, const vector<double>& w) {
    int n = v.size();
    vector<double> ret;
    ret.reserve(n);
    for (unsigned long i = 0; i != n; ++i) {
        ret.push_back(v[i] / w[i]);
    }
    return ret;
}



