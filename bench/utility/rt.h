#include<iostream>
#include<map>
#include<set>
#include<vector>
#include<fstream>
#include<random>
#include <bitset>

#include"rt_conf.h"

using namespace std;

#define DUMP_SYN_TRACE

//#define CLS 32
//#define DS 8
//#define TRACE
//#define PROFILE
//#define PROGILE_IBOUND

map<uint64_t, map<uint64_t, double>* > ref_mc;
uint64_t last_ref_id = -1;

vector< vector<double> > bit_mc(64, vector<double> (4, 0.0) );
uint64_t last_addr = -1;

map<uint64_t, vector< vector<double> >* > bit_mc_per_ref;
map<uint64_t, uint64_t> last_addr_per_ref;

uint64_t refT = 0;
map<uint64_t, uint64_t> lat;
map<uint64_t, double> ri;

void rtTmpAccess(uint64_t addr, uint64_t ref_id, uint64_t array_id, vector<int> idx) {
    uint64_t index = addr;
	addr = addr * DS / CLS;
    
    if (lat.find(addr) != lat.end()) {
        uint64_t singleRi = refT - lat[addr];
        if (ri.find(singleRi) == ri.end()) {
            ri[singleRi] = 1;
        } else {
            ri[singleRi] += 1;
        }
    }
    lat[addr] = refT;
    refT++;

    if (last_ref_id == -1) {
        last_ref_id = ref_id;
        last_addr = addr;
        last_addr_per_ref[ref_id] = addr;
        return;
    }

    if (ref_mc.find(last_ref_id) != ref_mc.end()) {
        if ((*ref_mc[last_ref_id]).find(ref_id) != (*ref_mc[last_ref_id]).end() ) {
            (*ref_mc[last_ref_id])[ref_id] ++;
        } else {
            (*ref_mc[last_ref_id])[ref_id] = 1;
        }
    } else {
        ref_mc[last_ref_id] = new map<uint64_t, double>;
        (*ref_mc[last_ref_id])[ref_id] = 1;
    }
 
    uint64_t bit_pos = 1;
    for (int i = 0; i < 64; i++) {
        uint64_t a = bit_pos & addr;
        uint64_t b = bit_pos & last_addr;
        if (a == 0 && b == 0) {
            bit_mc[i][0] += 1;       
        }
        if (a == 0 && b != 0) {
            bit_mc[i][1] += 1;
        }
        if (a != 0 && b == 0) {
            bit_mc[i][2] += 1;
        }
        if (a != 0 && b != 0) {
            bit_mc[i][3] += 1;
        }
        bit_pos = bit_pos << 1;
    }

    if (last_addr_per_ref.find(ref_id) != last_addr_per_ref.end() ) {
        if (bit_mc_per_ref.find(ref_id) == bit_mc_per_ref.end()) {
            bit_mc_per_ref[ref_id] = new vector< vector<double> >(64, vector<double> (4, 0.0) );
        }
        uint64_t bit_pos = 1;
        for (int i = 0; i < 64; i++) {
            uint64_t a = bit_pos & addr;
            uint64_t b = bit_pos & last_addr_per_ref[ref_id];
            if (a == 0 && b == 0) {
                (*bit_mc_per_ref[ref_id])[i][0] += 1;       
            }
            if (a == 0 && b != 0) {
                (*bit_mc_per_ref[ref_id])[i][1] += 1;
            }
            if (a != 0 && b == 0) {
                (*bit_mc_per_ref[ref_id])[i][2] += 1;
            }
            if (a != 0 && b != 0) {
                (*bit_mc_per_ref[ref_id])[i][3] += 1;
            }
            bit_pos = bit_pos << 1;
        }
    }

    last_ref_id = ref_id;
    last_addr = addr;
    last_addr_per_ref[ref_id] = addr;

    return;
}

uint64_t genAddr(uint64_t addr, vector< vector<double> >* bmc_used) {
    for (int i = 0; i < 64; i++) {
        uint64_t cur_bit = (1 << i) & addr;
        if (cur_bit == 0) {
            double flip_prob = double((*bmc_used)[i][1]) / ((*bmc_used)[i][0] + (*bmc_used)[i][1]);
            double rand_num = (double) rand() / (RAND_MAX);
            if (rand_num <= flip_prob) {
                addr |= (1 << i);
            }
        } else {
            double flip_prob = double((*bmc_used)[i][2]) / ((*bmc_used)[i][2] + (*bmc_used)[i][3]);
            double rand_num = (double) rand() / (RAND_MAX);
            if (rand_num <= flip_prob) {
                addr -= (1 << i);
            }
        }
    }
    return addr;
}

uint64_t filterAddrWithBmc(uint64_t addr, vector< vector<double> >* bmc_used) {
    uint64_t index = 0;
    for (auto elm : *bmc_used) {
        uint64_t mask = 1;
        mask = ((uint64_t) 1) << index;
        if (elm[0] == 1) {
            if ((addr & mask) != 0) { 
                addr -= mask;
            }
        }
        if (elm[4] == 1) {
            if ((addr & mask) == 0) { 
                addr |= mask;
            }
        }
        index += 1;
        //if (elm[0] == 1 || elm[4] == 1) cout << std::bitset<64>(addr) << " Index " << index << " addr " << endl;
    }
    return addr;
}

void predictionWithBmc(string conf = "") {
    map<uint64_t, uint64_t> lat_syn;
    map<uint64_t, double> ri_syn;
    uint64_t addr = ((uint64_t) rand() << 32) | rand();
    addr = filterAddrWithBmc(addr,  &bit_mc);    

#ifdef DUMP_SYN_TRACE
    cout << "SYN TRACE: ";
#endif
    for (int t = 0; t < 10000; t++) {
#ifdef DUMP_SYN_TRACE
        cout << addr << " ";
#endif
        if (lat_syn.find(addr) != lat_syn.end()) {
            uint64_t single_ri = t - lat_syn[addr];
            if (ri_syn.find(single_ri) != ri_syn.end()) {
                ri_syn[single_ri] += 1;
            } else {
                ri_syn[single_ri] = 1;
            }
        }
        lat_syn[addr] = t;
        addr = genAddr(addr, &bit_mc);
    }
#ifdef DUMP_SYN_TRACE
    cout << endl;
#endif


    fstream bmc_prediction_file;
    bmc_prediction_file.open("./result/trace_syn_ri/"+conf+".syn_ri", ios::out);
    double total_cnt = 0;
    for (auto elm : ri_syn) {
        total_cnt += elm.second;
    }
    for (auto elm : ri_syn) {
        bmc_prediction_file << elm.first << " " << elm.second / total_cnt << endl;
    }
    bmc_prediction_file.close();
}

uint64_t nextRef(uint64_t cur_ref) {
    double rand_num = (double) rand() / (RAND_MAX);
    for (auto ref_p : (*ref_mc[cur_ref])) {    
        if (rand_num <= ref_p.second) {
            return ref_p.first;
        } else {
            rand_num -= ref_p.second;
        }
        cur_ref = ref_p.first;
    }
    return cur_ref;
}

void predictionWithRefBmc(string conf = "") {
    map<uint64_t, uint64_t> lat_syn;
    map<uint64_t, double> ri_syn;
    
    set<uint64_t> refs;
    uint64_t cur_ref = -1;
    for (auto ref_ref_cnt : ref_mc) {
        if (refs.find(ref_ref_cnt.first) == refs.end()) {
            refs.insert(ref_ref_cnt.first);
            cur_ref = std::min(cur_ref, ref_ref_cnt.first);
        }
        for (auto ref_cnt : *(ref_ref_cnt.second)) {
            if (refs.find(ref_cnt.first) == refs.end()) {
                refs.insert(ref_cnt.first);
                cur_ref = std::min(cur_ref, ref_cnt.first);
            }
        }
    }
    map<uint64_t, uint64_t> addr_per_ref;
    for (auto ref : refs) {
        addr_per_ref[ref] = ((uint64_t) rand() << 32) | rand();
        addr_per_ref[ref] = filterAddrWithBmc(addr_per_ref[ref], bit_mc_per_ref[ref]);
    }
    
    for (int t = 0; t < 10000; t++) {
        if (lat_syn.find(addr_per_ref[cur_ref]) != lat_syn.end()) {
            uint64_t single_ri = t - lat_syn[addr_per_ref[cur_ref]];
            if (ri_syn.find(single_ri) != ri_syn.end()) {
                ri_syn[single_ri] += 1;
            } else {
                ri_syn[single_ri] = 1;
            }
        }
        lat_syn[addr_per_ref[cur_ref]] = t;
        addr_per_ref[cur_ref] = genAddr(addr_per_ref[cur_ref], bit_mc_per_ref[cur_ref]);
        cur_ref = nextRef(cur_ref);
    }
    fstream ref_bmc_prediction_file;
    ref_bmc_prediction_file.open("./result/trace_syn_ri/"+conf+".syn_ri_ref", ios::out);
    double total_cnt = 0;
    for (auto elm : ri_syn) {
        total_cnt += elm.second;
    }
    for (auto elm : ri_syn) {
        ref_bmc_prediction_file << elm.first << " " << elm.second / total_cnt << endl;
    }
    ref_bmc_prediction_file.close();
}

void dumpRIHistogram(string conf = "") {

    fstream bmc_file;
    bmc_file.open("./result/bmc/"+conf+".bmc", ios::out);
    
    uint64_t index = 0;
    for (auto & elm : bit_mc) {
        double sum = elm[0] + elm[1] + elm[2] + elm[3];
        elm[0] = elm[0]/sum;
        elm[1] = elm[1]/sum;
        elm[2] = elm[2]/sum;
        elm[3] = elm[3]/sum;
        bmc_file << index << ": " << elm[0] << " " << elm[1] << " " << elm[2] << " " << elm[3] << endl;
        index++;
    }
    bmc_file.close();    

    fstream refmc_file;
    refmc_file.open("./result/ref_bmc/"+conf+".refmc", ios::out);
    for (auto ref_ref_cnt : ref_mc) {
        double total_cnt = 0;
        for (auto ref_cnt : *(ref_ref_cnt.second)) {
            total_cnt += ref_cnt.second;
        }
        for (auto &ref_cnt : *(ref_ref_cnt.second)) {
            ref_cnt.second = ref_cnt.second / total_cnt;
            refmc_file << ref_ref_cnt.first << " -> " << ref_cnt.first << " : " << ref_cnt.second << endl;
        }
    }
    refmc_file.close();

    fstream ref_bmc_file;
    ref_bmc_file.open("./result/ref_bmc/"+conf+".ref_bmc", ios::out);
    for (auto bmc : bit_mc_per_ref) {
        ref_bmc_file << "ref " << bmc.first << endl;
        uint64_t index = 0;
        for (auto & elm : *bmc.second) {
            double sum = elm[0] + elm[1] + elm[2] + elm[3];
            elm[0] = elm[0]/sum;
            elm[1] = elm[1]/sum;
            elm[2] = elm[2]/sum;
            elm[3] = elm[3]/sum;
            ref_bmc_file << index << ": " << elm[0] << " " << elm[1] << " " << elm[2] << " " << elm[3] << endl;
            index++;
        }
    }
    ref_bmc_file.close();

    fstream ri_file;
    ri_file.open("./result/trace_ri/"+conf+".ri", ios::out);
    double cnt = 0;
    for (auto elm : ri) {
        cnt += elm.second;
    }
    for (auto elm : ri) {
        ri_file << elm.first << " " << elm.second / cnt << endl;
    }
    ri_file.close();
    return;
}
