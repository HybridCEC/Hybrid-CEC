#include "sweep.hpp"
#include "bitset.h"
  
extern "C" {
    #include "aiger.h"
}
void Sweep::aiger_preprocess() {
    aiger_read_from_file(aig, file);
    maxvar = aig->maxvar;
    output = aig->outputs[0].lit;
    inputs = aig->num_inputs;
    info.gates = aig->num_ands;
    info.vars = aig->maxvar;
    is_sat = false;
    
    q = new int[maxvar + 1];
    is_input = new int[maxvar + 1];
    used = new int[maxvar + 1];
    model = new int[maxvar + 1];
    bucket = new int[maxvar + 1];
    topo_counter = new int[maxvar + 1];
    fanout = new int[maxvar + 1];
    par = new int[maxvar + 1];
    active = new int[maxvar + 1];
    table = new int[maxvar + 1];
    cell = new type_cell[maxvar + 1];
    topo_index = new int[maxvar + 1];
    var_stamp = new int[maxvar + 1];
    xors_id = new int[maxvar + 1];
    xors_bot = new int[maxvar + 1];

    gate.push();
    for(int i=0; i<aig->num_ands; i++) {
        auto g = aig->ands[i];
        int o  = aiger_lit2var(g.lhs);
        int i0 = aiger_lit2var(g.rhs0);
        int i1 = aiger_lit2var(g.rhs1);
        gate.push();
        int id = gate.size() - 1;
        gate[id].push_in(i0 * (g.rhs0 & 1 ? -1 : 1));
        gate[id].push_in(i1 * (g.rhs1 & 1 ? -1 : 1));
        gate[id].out = o;
        cell[o] = (type_cell) {0, id, 0};
    }
    // printf("c INPUTS: %d, OUTPUTS: %d VARS: %d\n", inputs, aig->num_outputs, output >> 1);
    input = new int[inputs];
    for(int i=0; i<inputs; i++) {
        input[i] = aig->inputs[i].lit >> 1; 
    }
    for (int i = 1; i <= maxvar; i++) used[i] = is_input[i] = xors_id[i] = var_stamp[i] = bucket[i] = 0;
    for (int i = 1; i < gate.size(); i++) {
        used[abs(gate[i][0])] = 1;
        used[abs(gate[i][1])] = 1;
    }
    int k = 0;
    for (int i = 0; i < inputs; i++)
        if (used[input[i]])
            input[k++] = input[i];
        // else printf("c unused %d\n", input[i]);
    inputs =  k;
    for (int i = 0; i < inputs; i++) is_input[input[i]] = 1; 
    // printf("c USED INPUTS: %d\n", k);
}


Sweep* sweep_init(int argv, char* args[]) {
    // if(argv != 4) {
    //     printf("usage: ./sweep <SIMU-ROUNDS> <EPCEC-LIMITS>\n");
    //     exit(0);
    // }

    Sweep *S = new Sweep();
    S->aig = aiger_init();
    S->file = fopen(args[1], "r");
    // S->rounds = atoi(args[2]);
    S->rounds=20;
    return S;
}

void Sweep::solve() {
    aiger_preprocess();
    
    identify();

    cal_topo_index();

    cal_size();
    if (inputs <= 5) {
        seq_miter();
        print_result_and_exit("c find answer by miter!", nullptr);
    }
    else if (inputs <= 5) {
        // exact_pcec();
        epcec_in.clear();
        epcec_out.clear();
        for (int i = 0; i < inputs; i++) epcec_in.push(input[i]);
        epcec_out.push(output >> 1);
        int res = epcec(1);
        assert(!is_sat);
        print_result_and_exit("c find answer by epcec!", nullptr);
    }
    else {
        pcec();
        seq_CEC();
        assert(!is_sat);
        print_result_and_exit("c unsat proved by kissat!", nullptr);
    }
}

void Sweep::print_result_and_exit(const char* reason, int* input_vector) {
    // printf("%s\n", reason);
    if(!is_sat) {
        printf("Networks are equivalent\n");
    } else {
        printf("Networks are not equivalent\n");
        for (int i = 1; i <= maxvar; i++) table[i] = -1;

        for(int i=0; i<inputs; i++) {
            table[input[i]] = input_vector[i];
        }

        //simulate();

        for(int i=0; i<aig->num_inputs; i++) {
            int var = aig->inputs[i].lit >> 1;
            if(is_input[var]) {
                if(table[var]) {
                    printf("%d ", var);
                } else {
                    printf("%d ", -var);
                }
            } else {
                printf("%d ", var);
            }
        }
        puts("");
    }
    fflush(stdout);
    exit(0);
}

int main(int argv, char* args[]) {
    srand(19260817);
    Sweep *S = sweep_init(argv, args);
    S->solve();
    return 0;
}
