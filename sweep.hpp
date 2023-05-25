#ifndef sweep_hpp_INCLUDED
#define sweep_hpp_INCLUDED
#include <chrono>
#include <bits/stdc++.h>
#include <cstring>
#include "vec.hpp"

struct kissat;

struct aiger;
enum gate_type {And, Xor, UNUSE};
const std::string gate_type[6] = {"and", "xor", "majority", "HA", "FA", "UNUSE"};
class Bitset;

struct type_gate {
    vec<int> in;
    int ins, out, type;
    type_gate() : ins(0), out(0), type(0) {}
    void push_in (int x) { in .push(x); ins++; }
    int       &operator [] (int index)       { return in[index]; }
    const int& operator [] (int index) const { return in[index]; }
    type_gate& operator = (const type_gate& other) {
        ins = other.ins;
        out = other.out;
        type = other.type;
        in.growTo(ins);
        for (int i = 0; i < ins; i++) in[i] = other.in[i];
        return *this;
    }
};

struct type_cell {
    int fanouts = 0;
    int gate = 0;
    int is_carrier = 0;
};

struct cal_info {
    int vars = 0, gates = 0;
    int pairs = 0;
    int pairs_ec = 0, pairs_sat = 0;
    int pairs_ec_ver = 0, pairs_sat_ver = 0;
    int pairs_act = 0;
    int pairs_fanin = 0;
    int pairs_stru = 0;
    int pairs_ver = 0;
    double times_ec = 0, times_sat = 0;
};

class Sweep {
public:
    FILE* file;
    int maxvar, output, inputs, rounds, stamp = 0, xors = 0;
    double max_score = 0, max_connect = 0;
    bool is_sat;
    int *input, *topo_index, *used, *model, *size, *q, *active, *is_input, *topo_counter, *table, *var_stamp, *fanout, *par, *xors_bot, *xors_id, *bucket, *task_sign;
    vec<type_gate> gate;
    type_cell *cell;
    vec<int> *inv_C, *inv_V;
    int this_epcec = 0, simp_sz = 0;

    vec<int> epcec_in, epcec_out;
    std::chrono::high_resolution_clock::time_point clk_st = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<int, std::pair<int, int>>> equal_pairs;
    vec<std::pair<int, int>> verified_pairs;
    std::map<int, int> M;
    std::map<std::pair<int, int>, int> M_pair;
    aiger *aig;
    kissat * ksolver;
    cal_info info;

    void pre_seq_CEC();    
    void show_struct();
    int  find_par(int x);
    void cal_xor();
    void cal_size();
    void cal_topo_index();
    void identify();
    void recalculate_outs();
    bool match_xor(int x);
    int  merge(int x, int y);
    void fix_var(int x);

    void cal_xor_and_pis(int x, int y, int &pis, double &score);
    int  sign(int x);
    void print_model();
    void print_result_and_exit(const char* reason, int* input_vector);
    void solve();
    void aiger_preprocess();
    bool pcec();
    bool exact_pcec();
    bool epcec(bool exact = 1);
    bool _simulate(Bitset** result, void* _pool);
    bool check_var_struct(int x, int y);
    bool seq_CEC();
    bool seq_check_var_equal(int x, int y);
    bool seq_miter();

    void generate_input(Bitset** result, const vec<int> &input_list, void* _pool);
    void simulate();
    bool set_input_vector(int output, int value);
    bool try_flip_const(int var, int const_val);
};


#endif