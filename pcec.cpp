#include "sweep.hpp"
#include "aiger.h"
#include "bitset.h"
#include <unordered_map>

class BitsetPool {
public:
    std::queue<Bitset*> freeList;
    std::unordered_set<Bitset*> allocatedSet;
    int w;
    int total_mem_mb = 0;
    void expand() {
        const int n = 200;
        for(int i=0; i<n; i++) {
            Bitset* p = new Bitset();
            p->allocate(w);
            freeList.push(p);
        }
        // printf("c bitset pool expand n:%d w:%d\n", n, w);
        total_mem_mb += n * (w / 1024 / 1024 / 8);
    }
    void init(int _w) {
        w = _w;
        expand();
    }
    ~BitsetPool() {
        clear();
        while(!freeList.empty()) {
            auto p = freeList.front();
            freeList.pop();
            p->free();
            delete p;
        }
    }
    Bitset* allocate() {
        if(freeList.empty()) expand();
        Bitset* res = freeList.front();
        allocatedSet.insert(res);
        freeList.pop();
        return res;
    }
    void release(Bitset* p) {
        allocatedSet.erase(p);
        freeList.push(p);
    }
    void clear() {
        for(auto& p : allocatedSet) {
            freeList.push(p);
        }
        allocatedSet.clear();
    }
} global_pool;

bool Sweep::_simulate(Bitset** result, void* _pool) {

    BitsetPool *pool = (BitsetPool *) _pool;

    int *fanouts = new int[maxvar + 1];
    for (int i = 1; i <= maxvar; i++) {
        used[i] = 0;
        fanouts[i] = inv_C[i].size();
    }
    for (int i = 1; i <= gate.size(); i++) topo_counter[i] = 0;
    
    std::queue<int> q;
    
    for(int i = 0; i < epcec_in.size(); i++) {
        assert(active[epcec_in[i]]);
        q.push(epcec_in[i]);
        used[epcec_in[i]] = 2;
    }
    
    int o = output >> 1;

    while(!q.empty()) {
        int u = q.front();
        q.pop();

        int flag = 1;
        for (int i = 0; i < epcec_out.size(); i++)
            if (!used[epcec_out[i]]) {flag = 0; break;}
        if (flag) break;

        for(int i = 0; i < inv_C[u].size(); i++) {
            int c = inv_C[u][i];
            if (++topo_counter[c] != gate[c].ins) continue;
            int v = abs(gate[c].out);
            if (!active[v]) continue;
            if (var_stamp[v] != stamp) continue;
            
            q.push(v), used[v] = 1;
            result[v] = pool->allocate();
            if (gate[c].ins == 1) {
                int l1 = gate[c][0], v1 = abs(l1);
                result[v]->eqs(*result[v1], sign(gate[c].out) * sign(l1));
                fanouts[v1]--;
                if(fanouts[v1] == 0 && used[v1] != 2) pool->release(result[v1]);
                continue;    
            }
            int l1 = gate[c][0], l2 = gate[c][1];
            int v1 = abs(l1), v2 = abs(l2);
            if (gate[c].type == And)
                result[v]->ands(*result[v1], *result[v2], gate[c].out, l1, l2);
            else if (gate[c].type == Xor)
                result[v]->xors(*result[v1], *result[v2], gate[c].out, l1, l2);
            fanouts[v1]--;
            fanouts[v2]--;
            if(fanouts[v1] == 0 && used[v1] != 2) pool->release(result[v1]);
            if(fanouts[v2] == 0 && used[v2] != 2) pool->release(result[v2]);
        }
    }

    bool res = true;
    if (epcec_out.size() == 2) {
        Bitset* bitx = result[epcec_out[0]];
        Bitset* bity = result[epcec_out[1]];
        for (int i = 0; i < bitx -> m_size; i++)
            if (bitx->array[i] != bity->array[i]) {res = false; break;}
    } 
    else {
        Bitset& bit = *result[o];
        if (output & 1) bit.flip();
        for(int i = 0; i < bit.m_size; i++)
            if (bit.array[i] != 0) res = false;
        if(!res) {
            for (int i = 0; i < bit.m_size * 64; i++) {
                if (bit[i] == 0) continue;
                int input_vector[inputs];
                for(int j=0; j<inputs; j++) {
                    input_vector[j] = (*result[input[j]])[i];
                }
                is_sat = true;
                print_result_and_exit("c find answer by epcec !", input_vector);
            }
        }
    }
    delete fanouts;
    return res;
}

bool Sweep::epcec(bool exact) {

    // printf("c %cpcec checking with PI: ", exact ? 'e' : ' ');
    
    Bitset** result = new Bitset*[maxvar + 1];

    int num_inputs = epcec_in.size();

    const int maxR = 20;

    BitsetPool pool;
    pool.init(1LL << std::max(6, std::min(maxR, num_inputs)));

    if(num_inputs > maxR) {
        
        int extra_len = num_inputs - maxR;

        // 创建一个二进制数，用来枚举所有超过 maxR 位的，那些位的取值
        ull extra_values = 0;

        while(extra_values < (1LL << (extra_len)) ) {

            // printf("c epcec round [%d / %d]\n", extra_values, (1LL << (extra_len)));

            for(int i=0; i<num_inputs; i++) {
                result[epcec_in[i]] = pool.allocate();
            }

            for(int i=0; i<extra_len; i++) {
                int input_var = epcec_in[i];
                ull val = extra_values & (1LL << (extra_len - i - 1));
                if(val != 0) val = ~(val = 0);

                for(int j=0; j<result[input_var]->m_size; j++) {
                    result[input_var]->array[j] = val;
                }
            }
            int sz = 1 << maxR;
            int unit = sz;
            for(int i=extra_len; i<num_inputs; i++) {
                int input_var = epcec_in[i];
                unit >>= 1;
                if(unit >= 64) {
                    const ull all_zero = 0;
                    const ull all_one = ~all_zero;
                    for(int j=0; j<sz/64; j++) {
                        int value = (j * 64 / unit) % 2;
                        assert(j < result[input_var]->m_size);
                        if(value) result[input_var]->array[j] = all_one;
                        else result[input_var]->array[j] = all_zero;
                    }
                } else {
                    for(int j=0; j<64; j++) {
                        int value = (j >> (num_inputs - i - 1)) & 1;
                        if(value) result[input_var]->reset(j);
                        else result[input_var]->set(j);
                    }
                    
                    for(int j=1; j<sz/64; j++) {
                        result[input_var]->array[j] = result[input_var]->array[0];
                    }
                }
            }

            int res = _simulate(result, &pool);
            pool.clear();

            extra_values++;
            if (!res) return 0;
        }
    } else if(num_inputs > 6) {
        for(int i=0; i<num_inputs; i++) {
            result[epcec_in[i]] = pool.allocate();
        }
        int sz = 1 << epcec_in.size();
        int unit = sz;
        for(int i=0; i<epcec_in.size(); i++) {
            int input_var = epcec_in[i];
            unit >>= 1;
            if(unit >= 64) {
                const ull all_zero = 0;
                const ull all_one = ~all_zero;
                for(int j=0; j<sz/64; j++) {
                    int value = (j * 64 / unit) % 2;
                    if(value) result[input_var]->array[j] = all_one;
                    else result[input_var]->array[j] = all_zero;
                }
            } else {
                for(int j=0; j<64; j++) {
                    int value = (j >> (num_inputs - i - 1)) & 1;
                    if(value) result[input_var]->reset(j);
                    else result[input_var]->set(j);
                }
                
                for(int j=1; j<sz/64; j++) {
                    result[input_var]->array[j] = result[input_var]->array[0];
                }
            }
        }
        return _simulate(result, &pool);
    } else {
        assert(num_inputs <= 6);
        for(int i=0; i<num_inputs; i++) {
            result[epcec_in[i]] = pool.allocate();
        }


        // printf("\n");
        for(int i=0; i<num_inputs; i++) {
            assert(result[epcec_in[i]]->m_size == 1);

            ull s = 0;

            for(int j=0; j<(1LL<<num_inputs); j++) {
                s <<= 1;
                if(j & (1 << (num_inputs-i-1))) {
                    s |= 1;
                }
            }
            result[epcec_in[i]]->array[0] = s;
        }

        // printf("\n");

        // for(int i=0; i<num_inputs; i++) {
        //      result[epcec_in[i]]->print();
        // }
        // exit(0);

        return _simulate(result, &pool);
    }
    return true;

}

bool Sweep::pcec() {
    auto clk_st = std::chrono::high_resolution_clock::now();
    int R, w;
    if(rounds > 20) {
        rounds = 1 << (rounds - 20);
        R = 20;
    } else {
        R = rounds;
        rounds = 1;
    }
    w = 1ll << R;

    global_pool.init(w);
    Bitset** result = new Bitset*[maxvar + 1];
    
    vec<int> zero, one;
    std::queue<int> q;
    
    Bitset A, B;  //A = 0, B = 1;
    A.allocate(w); A.set(); A.hash();
    B.allocate(w); B.reset(); B.hash();
    ull one_val = A.hashval;
    ull zero_val = B.hashval;
    // printf("%llu %llu\n", zero_val, one_val);
    int o = output >> 1;

    std::vector<std::vector<int>> classes;
    classes.push_back(std::vector<int>());

    for(int i=1; i<=maxvar; i++) {
        int c = cell[i].gate;
        if (!c || gate[c].type == UNUSE) continue;
        classes[0].push_back(i);
    }

    int *fanouts = new int[maxvar + 1];
    
    for(int rd=0; rd<rounds; rd++) {
        std::unordered_map<int, ull> HashV;
        
        for (int i = 1; i <= maxvar; i++) {
            topo_counter[i] = 0;
            fanouts[i] = inv_V[i].size();
        }

        global_pool.clear();

        // printf("c random test round [%d / %d]\n", rd + 1, rounds);
        for(int i=0; i<inputs; i++) {
            int input_var = input[i];
            result[input_var] = global_pool.allocate();
            result[input_var]->random();
            result[input_var]->hash();
            HashV[input_var] = result[input_var]->hashval;
            q.push(input_var);
        }
        // printf("c time1: %d\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - clk_st).count());
        while(!q.empty()) {     
            int u = q.front();
            q.pop();
            for(int i = 0; i < inv_C[u].size(); i++) {
                int c = inv_C[u][i];
                int v = abs(gate[c].out);
                topo_counter[c]++;
                assert(gate[c].ins == 2);
                if(topo_counter[c] == gate[c].ins) {
                    q.push(v);
                    result[v] = global_pool.allocate();
                    int l1 = gate[c][0], l2 = gate[c][1];
                    int v1 = abs(l1), v2 = abs(l2);

                    // try to release mem space
                    fanouts[v1]--;
                    fanouts[v2]--;
                    if(fanouts[v1] == 0 && !is_input[v1]) global_pool.release(result[v1]);
                    if(fanouts[v2] == 0 && !is_input[v2]) global_pool.release(result[v2]);

                    if (gate[c].type ==  And)
                        result[v]->ands(*result[v1], *result[v2], gate[c].out, l1, l2);
                    else if (gate[c].type ==  Xor)
                        result[v]->xors(*result[v1], *result[v2], gate[c].out, l1, l2);
                    result[v]->hash();
                    ull hval = result[v]->hashval;
                    HashV[v] = hval;
                    if (hval == zero_val && v != o) {
                        if (*result[v] == B) zero.push(v);
                    }
                    else if (hval == one_val && v != o) {
                        if (*result[v] == A) one.push(v);
                    }
                }
                
            }
        }

        // printf("c classes-size: %d\n", classes.size());

        std::vector<std::vector<int>> _classes;
        for(int i=0; i<classes.size(); i++) {
            std::vector<int> &c = classes[i];
            std::unordered_map<ull, std::vector<int>> mp;
            for(int j=0; j<c.size(); j++) {
                int var = c[j];
                mp[HashV[var]].push_back(var);
            }

            for(auto& [k, v] : mp) {
                if(v.size() <= 1) continue;
                if(v.size() >= 20) {
                    // printf("c skip too large class!");
                    continue;
                }
                _classes.push_back(v);
            }
        }

        classes = _classes;

        Bitset& res = *result[o];
        if (output & 1) res.flip();
        for(int i=0; i<w; i++) {
            if(res[i] == 0) continue;
            int input_vector[inputs];
            for(int j=0; j<inputs; j++) {
                input_vector[j] = (*result[input[j]])[i];
            }
            is_sat = true;
            print_result_and_exit("c find answer by pcec!", input_vector);
        }
    }

    analyze_result:

    // printf("val = 0: ");
    // for (int i = 0; i < zero.size(); i++) printf("%d ", zero[i]);
    // printf("\nval = 1: ");
    // for (int i = 0; i < one.size(); i++) printf("%d ", one[i]);
    // puts("");

    // for(int i=0; i<zero.size(); i++) try_flip_const(zero[i], 0);
    // for(int i=0; i<one.size(); i++) try_flip_const(one[i], 1);

    // printf("c time2: %d\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - clk_st).count());
    
    // printf("c bitset pool mem: %d MB\n", global_pool.total_mem_mb);
    
    delete result;
    delete fanouts;

    // printf("c equivalent classes: ");
    // for(auto& clas : classes) {
    //     printf("{ ");
    //     for(auto&v : clas) {
    //         printf("%d ", v);
    //     }
    //     printf("}\n");
    // }

    for(auto& V : classes) {
        for (int j = 0; j < V.size(); j++) {
            for(int k = j+1; k <V.size(); k++) {
                equal_pairs.push_back(std::make_pair(std::max(topo_index[V[j]], topo_index[V[k]]), std::make_pair(V[j], V[k])));
            }
        }
    }
    
    sort(equal_pairs.begin(), equal_pairs.end());
    // printf("c num of equivalent pairs: %d\n", equal_pairs.size());
    return true;
}

// void Sweep::pre_seq_CEC() {
//     double score;
//     int pis;
//     for (int i = equal_pairs.size() - 1; i >= 0; i--) {
//         int x = equal_pairs[i].second.first, y = equal_pairs[i].second.second;
//         if (!active[x] || !active[y]) {
//             task_sign[i] = 0;
//             printf("c [%d/%d] trying to prove equivalance of pair (%d, %d) (size: %d %d)\n", i+1, equal_pairs.size(), x, y, size[x], size[y]); 
//             printf("c can skip cs not active +++++++++++++++++++++++++++\n");
//             skip_act++; continue;
//         }
//         cal_xor_and_pis(x, y, pis, score);
//         if (pis <= 20) task_sign[i] = 0;
//         // if (pis <= 32 && score >= 0.32) task_sign[i] = 0;
//         //TODO: if fanin of x or y are in equal_pairs, set task_sign[i] to -1
//         if (task_sign[i] == -1) continue; 
//         printf("c [%d/%d] trying to prove equivalance of pair (%d, %d) (size: %d %d)\n", i+1, equal_pairs.size(), x, y, size[x], size[y]);     
//         printf("c inputs number: %d\n", pis); 
//         int mx = M[size[x]], my = M[size[y]];
//         int c1 = cell[x].gate, c2 = cell[y].gate, res;
//         if (mx == my && mx) {
//             int px = verified_pairs[mx].first;
//             int py = verified_pairs[mx].second;
//             if ((check_var_struct(x, px) && check_var_struct(y, py)) || 
//                 (check_var_struct(y, px) && check_var_struct(x, py))) {
//                 printf("c can skip cs struct equal =============================\n");   
//                 skip_struct++;
//                 goto equal;
//             }
//         }
//         epcec_out.clear();
//         epcec_out.push(abs(x));
//         epcec_out.push(abs(y));
//         res = epcec(1);
//         if (!res) continue;
// equal:
//         eqs++;
//         merge(x, y);
//         M[size[x]] = verified_pairs.size();
//         M[size[y]] = verified_pairs.size();
//         verified_pairs.push(std::make_pair(x, y));
//     }
//     //can simplify ?
// }

bool Sweep::seq_CEC() {    
    // printf("use time: %d\n\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - clk_st).count());    
    // printf("c sequential mode\n");
    // printf("c checkpair (%d, %d) (size: %d %d)\n", 992, 993, size[992], size[993]);
    // int resc = seq_check_var_equal(992, 993);
    // printf("res is %d, use time: %d\n\n", resc, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - clk_st).count());  
    info.pairs = equal_pairs.size();  
    task_sign = new int[equal_pairs.size()];
    for (int i = 0; i < equal_pairs.size(); i++) task_sign[i] = -1;
    // pre_seq_CEC();
    for(int i=0; i<equal_pairs.size(); i++) {
        auto clk_ste = std::chrono::high_resolution_clock::now();
        if (task_sign[i] != -1) continue;
        this_epcec = 0;
        auto& pr = equal_pairs[i];
        int eq = -1;
        int x = pr.second.first, y = pr.second.second;
        // printf("c [%d/%d] trying to prove equivalance of pair (%d, %d) (size: %d %d)\n", i+1, equal_pairs.size(), x, y, size[x], size[y]);
        if (!active[x] || !active[y]) {
            // printf("c can skip cs not active +++++++++++++++++++++++++++\n");
            info.pairs_act++;
            continue;
        }
        int mx = M[size[x]], my = M[size[y]];
        int c1 = cell[x].gate, c2 = cell[y].gate;
        if (mx == my && mx) {
            int px = verified_pairs[mx].first;
            int py = verified_pairs[mx].second;
            if ((check_var_struct(x, px) && check_var_struct(y, py)) || 
                (check_var_struct(y, px) && check_var_struct(x, py))) {
                info.pairs_stru++;
                eq = 1;
                // printf("c can skip cs struct equal =============================\n");
                goto verified;
            }
        }

        if (gate[c1].type == gate[c2].type && gate[c1].ins) {
            int a = gate[c1][0], b = gate[c1][1], c = gate[c2][0], d = gate[c2][1];
            if (a == c && b == d) eq = 1;
            if (a == d && b == c) eq = 1;

            if (eq == 1) {
                info.pairs_fanin++;
                // printf("c can skip cs fanin equal ******************************\n");
                goto verified;
            }

        }
        eq = seq_check_var_equal(x, y);
        // printf("c use SAT/epcec proving result(is equal): %d\t", eq);


verified:
        double use_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - clk_ste).count() * 0.001;
        if (this_epcec == 1) {
            info.pairs_ec++;
            info.times_ec += use_time;
        }
        else if (this_epcec == 2) {
            info.pairs_sat++;
            info.times_sat += use_time;
        }
        if(eq == 1) {
            info.pairs_ver++;
            if (this_epcec == 1) info.pairs_ec_ver++; else if (this_epcec == 2) info.pairs_sat_ver++;
            M[size[x]] = verified_pairs.size();
            M[size[y]] = verified_pairs.size();
            verified_pairs.push(std::make_pair(x, y));
            // M_pair[std::make_pair(x, y)] = 1;
            // M_pair[std::make_pair(y, x)] = 1;
            // assert(x > 0 && y > 0);
            simp_sz += merge(x, y);
        }
        // printf("simplified: %d\n\n", simp_sz);    
    }
    // printf("c max score %.2lf\n", max_score);
    // printf("c max connected %.2lf\n", max_connect);
    // printf("c vars (gates) %d %d\n", info.vars, info.gates);
    // printf("c pairs (non-active, struct-equal, same-fanin, epcec, sat-solver) %d %d %d %d %d %d\n", info.pairs, info.pairs_act, info.pairs_stru, info.pairs_fanin, info.pairs_ec, info.pairs_sat);
    // printf("c use time (epcec, sat-solver) %.2lf %.2lf\n", info.times_ec, info.times_sat);
    // printf("c verified pairs (epcec, sat-solver) %d %d %d\n", info.pairs_ver, info.pairs_ec_ver, info.pairs_sat_ver);
    
    int res = seq_miter();
    // printf("c entire-problem: %d\n", res);
    // printf("c miter time: %d\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - clk_st).count());

    return res;
}
