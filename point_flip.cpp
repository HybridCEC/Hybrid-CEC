#include "sweep.hpp"

extern "C" {
   #include "src/internal.h"
   #include "src/witness.h"
}

void Sweep::simulate() {

    std::queue<int> q;

    for(int i=0; i<inputs; i++) {
        q.push(input[i]);
    }

    for (int i = 1; i <= maxvar; i++) topo_counter[i] = 0;

    while(!q.empty()) {
        int u = q.front();
        q.pop();

        for(int i=0; i<inv_V[u].size(); i++) {
            int v = inv_V[u][i];
            topo_counter[v]++;
            if(topo_counter[v] == 2) {
                q.push(v);
                int c = cell[v].gate;
                printf("V: %d (%d) max:%d c:%d\n", v, gate[c].ins, maxvar, c);

                assert(gate[c].ins == 2);
                printf("in0: %d in1: %d\n", gate[c][0], gate[c][1]);
                assert(gate[c].type != UNUSE);

                int rhs0 = gate[c].in[0] < 0 ? !table[abs(gate[c].in[0])] : table[abs(gate[c].in[0])];
                int rhs1 = gate[c].in[1] < 0 ? !table[abs(gate[c].in[1])] : table[abs(gate[c].in[1])];

                if(gate[c].type == And) {
                    table[v] = rhs0 & rhs1;
                } else {
                    assert(gate[c].type == Xor);
                    table[v] = rhs0 ^ rhs1;
                }
                
            }
        }
    }
}

bool Sweep::set_input_vector(int output, int value) {
    
    ksolver = kissat_init();

    std::queue<int> q;
    for (int i = 1; i <= maxvar; i++) used[i] = 0;

    q.push(output);
    used[output] = true;

    while(!q.empty()) {
        int u = q.front();
        q.pop();
        auto& gt = gate[cell[u].gate];
        if(gt.ins == 0) continue;
        assert(gt.ins == 2);
        int a = gt[0], b = gt[1];
        if(gt.type == And) {
            kissat_add(ksolver, -u); kissat_add(ksolver, a); kissat_add(ksolver, 0);
            kissat_add(ksolver, -u); kissat_add(ksolver, b); kissat_add(ksolver, 0);
            kissat_add(ksolver, -a); kissat_add(ksolver, -b); kissat_add(ksolver, u); kissat_add(ksolver, 0);
        } else {
            assert(gt.type == Xor);
            kissat_add(ksolver, -a); kissat_add(ksolver, -b); kissat_add(ksolver, -u); kissat_add(ksolver, 0);
            kissat_add(ksolver, -a); kissat_add(ksolver, b); kissat_add(ksolver, u); kissat_add(ksolver, 0);
            kissat_add(ksolver, a); kissat_add(ksolver, -b); kissat_add(ksolver, u); kissat_add(ksolver, 0);
            kissat_add(ksolver, a); kissat_add(ksolver, b); kissat_add(ksolver, -u); kissat_add(ksolver, 0);
        }
        
        a = abs(a), b = abs(b);
        if(!used[a]) q.push(a), used[a] = true;
        if(!used[b]) q.push(b), used[b] = true;
    }

    if(value) {
        kissat_add(ksolver, output); kissat_add(ksolver, 0);
    } else {
        kissat_add(ksolver, -output); kissat_add(ksolver, 0);
    }

    int result = kissat_solve(ksolver);

    if(result == 20) {
        delete ksolver;
        return false;
    }
 
    for(int i=0; i<inputs; i++) {
        int var = input[i];
        int tmp = kissat_value(ksolver, var);
        if (!tmp) tmp = var;
        if(tmp > 0) table[var] = 1;
        else table[var] = 0;
    }

    delete ksolver;
    return true;
}

bool Sweep::try_flip_const(int var, int const_val) {

    printf("c VAR-FLIP: try to flip const var %d = %d\n", var, const_val);

    for (int i = 1; i <= maxvar; i++) table[i] = -1;
    

    // get counter-example
    bool success = set_input_vector(var, !const_val);

    // is const var
    if(!success) {
        printf("c find const var: %d = %d\n", var, const_val);
        return true;
    }

    simulate();

    assert(table[var] == !const_val);

    int o = output >> 1;

    if(table[o] != (output & 1)) {
        int input_vector[inputs];
        for(int i=0; i<inputs; i++) {
            input_vector[i] = table[input[i]];
        }
        is_sat = true;
        print_result_and_exit("c find answer by var-flipping", input_vector);
    }

    return false;
}