/* 
 * File:   ODD_delayed_values.h
 * Author: luecken
 *
 * Created on 3. März 2014, 14:19
 */

#ifndef ODD_DELAYED_VALUES_H
#define	ODD_DELAYED_VALUES_H

#include <map>
#include <vector>
#include <memory>
#include <iostream>
#include <list>
#include <algorithm>
#include "ODD_utils.h"

/*
 * ODD_delayed_values provides a data structure which 
 * has two index sets: one is the original indexing of delays provided
 * by the user and the other is the indexing obtained after merging equal delays
 * and sorting the delays (descendingly);
 * the indices actually used can be switched by useUserIndices(bool val);
 * they should be switched on each time the user supplied rhs function is called,
 * since in its scope, they will be valid;
 * for internal processing (e.g. interpolation), the optimized indices should be used.
 * 
 */

class ODD_delayed_values {
    std::vector<std::map<size_t, double> > data; // underlying data
    std::vector<size_t> ui; // mapping from user indices to sorted indices 
    bool use_user_indices = false;
    // indicates from which index i on the delays are 0
    size_t tau0_offset;
public:
    typedef std::vector<std::map<size_t, double> >::iterator iterator;
    typedef std::vector<std::map<size_t, double> >::const_iterator const_iterator;

    ODD_delayed_values() {
    };

    // make the user index mapping
    //    void init(std::shared_ptr<const std::vector<size_t> > mapping, size_t tau0_offset) {

    void init(const std::vector<size_t>* mapping, size_t tau0_offset) {
        ui = *mapping;
        size_t data_size = 0;
        if (mapping->size() != 0) data_size = *std::max_element(std::begin(ui), std::end(ui)) + 1;
        std::cout << "ODD_delayed_values::data size = " << data_size << std::endl;
        data = std::vector<std::map<size_t, double> >(data_size);
        this->tau0_offset = tau0_offset;
    }

    std::map<size_t, double>& operator[] (size_t index) {
        if (use_user_indices) return data[ui[index]];
        else return data[index];
    }

    size_t getTau0_offset() const {
        return tau0_offset;
    }

    // fill delays tau_i==0 with the value from the provided current state vector x_current

    void fillZeroDelays(const std::vector<double>& x_current) {
        for (size_t j = tau0_offset; j != data.size(); ++j) {
            for (auto& e : data[j]) {
                e.second = x_current[e.first];
            }
        }
    }

    void useUserIndices(bool val) {
        use_user_indices = val;
    }

    size_t size() const {
        if (use_user_indices) return ui.size();
        else return data.size();
    }

    void clear() {
        for (auto e : data) e.clear();
    }

    iterator begin() {
        return data.begin();
    }

    iterator end() {
        return data.end();
    }

    const_iterator begin() const {
        return data.begin();
    }

    const_iterator end() const {
        return data.end();
    }

    void print() const {
        std::cout << "delayed values: [";
        for (auto e : data) {
            for(auto p : e) {
            std::cout << p.second << ", ";
        }
        }
        std::cout << "]" << std::endl;
    }
};

/* create and return a mapping from user indices to a set of reduced indices
 *      p: {0,..., max(keys(m))-1} -> {0,..., set(values(m)).size()-1} 
 * that serves as a view from the original indices in m to a decreasingly sorted vector w
 * which contains each element once, i.e. w[p[j]]=m[j]
 */
template<typename T>
void sort_and_make_index_map(const std::map<size_t, T>& m, std::vector<T>& sort_target, std::vector<size_t>& imap) {
    std::cout << "sort_and_make_index_map()" << std::endl;
    //    std::cout << "size(m) = " << m.size() << std::endl;
    auto delay_indices = getKeys(m);
    auto max_index = *std::max_element(begin(delay_indices), end(delay_indices));
    //    std::cout << "maximal index in m = " << max_index << std::endl;

    imap = std::vector<size_t > (max_index + 1);
    size_t Nd = nr_of_distinct_elements(m);
    sort_target.clear();
    sort_target.reserve(Nd);

    // make lists that contain all indices with same elements
    // list l[t] will contain all keys from m with value t. Note that l will be sorted ascendingly
    std::map<T, std::list<size_t> > l;
    for (auto e : m) {
        //        std::cout << "t_" << e.first << " = " << e.second << std::endl;
        l[e.second].push_back(e.first);
    }

    // from l make the index mapping and the sorted value list
    for (auto e = l.rbegin(); e != l.rend(); ++e) {
        sort_target.push_back(e->first);
        size_t j = sort_target.size() - 1;
        for (auto index : e->second) {
            imap[index] = j;
        }
    }

    std::cout << "imap:" << std::endl;
    for (size_t j = 0; j != imap.size(); ++j) {
        std::cout << j << " ";
    }
    std::cout << "\n";
    for (size_t j = 0; j != imap.size(); ++j) {
        std::cout << imap[j] << " ";
    }
    std::cout << std::endl;

}


#endif	/* ODD_DELAYED_VALUES_H */

