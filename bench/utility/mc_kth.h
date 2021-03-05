#include<iostream>
#include<map>
#include<set>
#include<vector>
#include<fstream>
#include<random>
#include<bitset>

#include"rt_conf.h"

using namespace std;

//#define DUMP_KTM_CONVERSION
//#define DUMP_KTM_TO_BIT_MAPPING
//#define DUMP_BMC
//#define DUMP_SYN_TRACE

vector< map<uint64_t, pair<double, double> >* > bit_mc(64, nullptr);
vector< uint64_t > last_k_addrs;

uint64_t refT = 0;

map<uint64_t, uint64_t> lat;
map<uint64_t, double> ri;

// This covert KLEVEL * TLEVEL bits to an uint64_t
uint64_t covertKTBitMatrixToInt(vector< uint64_t >* pre_k_addrs, int ith_bit) {
    uint64_t converted_value = 0;
    int curr_index = 0;
#ifdef DUMP_KTM_CONVERSION
    if (ith_bit == 0 || ith_bit == 1) {
        cout << "Pre K addrs: ";
        for (auto elm : *pre_k_addrs) { cout << bitset<64>(elm).to_string() << " "; }
    }
#endif
    for (int i = 0; i < KLEVEL; i++) {
        for (int j = 0; j < TLEVEL; j++) {
            
            bool curr_bit = ((((uint64_t) 1) << ((ith_bit - i) % 64))) & (*pre_k_addrs)[j];
#ifdef DUMP_KTM_CONVERSION
            if (ith_bit == 0 || ith_bit == 1) { cout << "B" << ((ith_bit - i) % 64) << "T" << j << " " << curr_bit << " "; }
#endif
            if (curr_bit) {
                converted_value += ((uint64_t) 1) << curr_index;
            }
            curr_index += 1;
        }
    }
#ifdef DUMP_KTM_CONVERSION
    if (ith_bit == 0 || ith_bit == 1) { cout << "Converted: " << converted_value << endl; }
#endif

    return converted_value;
}

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
    
    if (last_k_addrs.size() < TLEVEL) {
        last_k_addrs.push_back(addr);
        return;
    }
    
#ifdef DUMP_KTM_TO_BIT_MAPPING
    cout << "K addrs : ";
    for (auto elm : last_k_addrs) { cout
        << bitset<64>(elm).to_string() << " "; }
    cout << endl;
#endif
    
    for (int i = 0; i < 64; i++) {
        uint64_t encodedKTBMatrix = covertKTBitMatrixToInt( &last_k_addrs, i );
        
        bool cur_bit = addr & (((uint64_t) 1) << i);

#ifdef DUMP_KTM_TO_BIT_MAPPING
        cout << "B" << i << " " << encodedKTBMatrix << " " << cur_bit << " ";
#endif
        
        if (bit_mc[i] == nullptr) {
            bit_mc[i] = new map<uint64_t, pair<double, double> >;
        }
        
        if (bit_mc[i]->find(encodedKTBMatrix) == bit_mc[i]->end()) {
            (*bit_mc[i])[encodedKTBMatrix] = make_pair(0, 0);
        }
        
        if (cur_bit) {
            (*bit_mc[i])[encodedKTBMatrix].second += 1;
        } else {
            (*bit_mc[i])[encodedKTBMatrix].first += 1;
        }
    }
    
#ifdef DUMP_KTM_TO_BIT_MAPPING
    cout << endl;
#endif
    
    last_k_addrs.erase(last_k_addrs.begin());
    last_k_addrs.push_back(addr);
    
    return;
}

uint64_t genAddr(vector<uint64_t>* pre_k_addrs, vector< map<uint64_t, pair<double, double> >* >* bmc_used) {
    uint64_t addr_new = 0;
    for (int i = 0; i < 64; i++) {
        
        uint64_t encodedKTBMatrix = covertKTBitMatrixToInt( pre_k_addrs , i);
        
        if ((*bmc_used)[i]->find(encodedKTBMatrix) != (*bmc_used)[i]->end()) {
            double rand_num = (double) rand() / (RAND_MAX);
            if (rand_num > (*(*bmc_used)[i])[encodedKTBMatrix].first) {
                addr_new = addr_new | (((uint64_t) 1) << i);
            }
        } else {
            auto it = (*bmc_used)[i]->begin();
            std::advance(it, rand() % (*bmc_used)[i]->size());
            uint64_t random_key = it->first;
            double rand_num = (double) rand() / (RAND_MAX);
            if (rand_num > (*(*bmc_used)[i])[random_key].first) {
                addr_new = addr_new | (((uint64_t) 1) << i);
            }
        }
    }
    return addr_new;
}

void normalizeBMC(vector< map<uint64_t, pair<double, double> >* >* bmc_used) {
    for (auto map_ptr : *bmc_used) {
        for (auto elm : *map_ptr) {
            double sum = elm.second.first + elm.second.second;
#ifdef DUMP_BMC
            cout << elm.first << " " << (*map_ptr)[elm.first].first << " " << (*map_ptr)[elm.first].second << endl;
#endif
            (*map_ptr)[elm.first] = make_pair(((double) elm.second.first) / sum, ((double)elm.second.second)/ sum);
#ifdef DUMP_BMC
            cout << "Normalized: " << elm.first << " " << (*map_ptr)[elm.first].first << " " << (*map_ptr)[elm.first].second << endl;
#endif
        }
    }
}

void predictionWithBmc(string conf = "") {
    map<uint64_t, uint64_t> lat_syn;
    map<uint64_t, double> ri_syn;
    
    //normalizeBMC(&bit_mc);
    
    vector<uint64_t> pre_t_addrs;
    for (int i = 0; i < TLEVEL; i++) {
        pre_t_addrs.push_back(((uint64_t) rand() << 32) | rand());
    }
    
    for (auto elm : pre_t_addrs) {
        cout << elm << " ";
    }
    cout << endl;
    
    uint64_t addr = pre_t_addrs.back();
    
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
        addr = genAddr(&pre_t_addrs, &bit_mc);
        pre_t_addrs.erase(pre_t_addrs.begin());
        pre_t_addrs.push_back(addr);
    }
#ifdef DUMP_SYN_TRACE
    cout << endl;
#endif

    fstream bmc_prediction_file;
    bmc_prediction_file.open("./result/trace_syn_ri/"+conf+".syn_ri_" + to_string(KLEVEL) + "_" + to_string(TLEVEL), ios::out);
    double total_cnt = 0;
    for (auto elm : ri_syn) {
        total_cnt += elm.second;
    }
    for (auto elm : ri_syn) {
        bmc_prediction_file << elm.first << " " << elm.second / total_cnt << endl;
    }
    bmc_prediction_file.close();
}

void dumpRIHistogram(string conf = "") {

    fstream bmc_file;
    bmc_file.open("./result/bmc/"+conf+".bmc_" + to_string(KLEVEL) + "_" + to_string(TLEVEL), ios::out);
    
    normalizeBMC(&bit_mc);
    
    uint64_t index = 0;
    for (auto map_ptr : bit_mc) {
        bmc_file << index << ": " << endl;
        for (auto elm : *map_ptr) {
            bmc_file << "    "  << elm.first << " " << elm.second.first << " " << elm.second.second << endl;
        }
        index++;
    }
    bmc_file.close();

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
