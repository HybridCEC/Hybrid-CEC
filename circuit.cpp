#include "sweep.hpp"

int Sweep::sign(int x) {
    return x ? (x > 0 ? 1 : -1) : 0;
}

void Sweep::cal_xor_and_pis(int x, int y, int &pis, double &score) {
    pis = 0, score = 0;
    int l = 1, r = 0;
    q[++r] = abs(x), used[abs(x)] = true;
    q[++r] = abs(y), used[abs(y)] = true;
    var_stamp[abs(x)] = var_stamp[abs(y)] = ++stamp;
    epcec_in.clear();
    while(l <= r) {
        int u = q[l++];
        if (is_input[u]) epcec_in.push(u), ++pis;
        auto& gt = gate[cell[u].gate];
        if(gt.ins == 0) continue;
        assert(gt.ins == 2);
        //need fixed
        int a = abs(gt[0]), b = abs(gt[1]);
        assert(active[a] && active[b]);
        if (var_stamp[a] != stamp) q[++r] = a, var_stamp[a] = stamp;
        if (var_stamp[b] != stamp) q[++r] = b, var_stamp[b] = stamp;
    }
    for (int i = 1; i <= r; i++) par[q[i]] = q[i], xors_bot[q[i]] = 0;
    for (int i = 1; i <= r; i++) {
        int u = q[i];
        if (gate[cell[u].gate].type != Xor) continue;
        //need fixed ins = 2
        int a = abs(gate[cell[u].gate][0]), b = abs(gate[cell[u].gate][1]);
        int c1 = cell[a].gate, c2 = cell[b].gate;
        if (!c1 || !c2 || gate[c1].type != Xor || gate[c2].type != Xor) continue;
        par[find_par(a)] = find_par(u);
        par[find_par(b)] = find_par(u);
    }
    for (int i = 1; i <= r; i++) {
        int u = q[i];
        if (find_par(u) == u) continue;
        xors_bot[par[u]]++;
    }
    // printf("c xors block size: ");
    int mx = 0;
    for (int i = 1; i <= r; i++) {
        int u = q[i];
        if (xors_bot[u] <= 2) continue;
        bucket[++xors_bot[u]]++;
        mx = std::max(mx, xors_bot[u]);
        // printf("%d ", xors_bot[u]);
    }
    // puts("");
    score = mx;
    for (int i = 3; i < mx; i++) {
        bucket[i + 1] += bucket[i] >> 2;
        bucket[i] %= 2;
    }
    while (bucket[mx] >= 4) score++, bucket[mx] /= 4;
    for (int i = 1; i <= mx; i++) bucket[i] = 0;
    score /= 1.0 * pis;
    if (pis > 20) max_score = std::max(max_score, score);
}

void Sweep::cal_size() {
    size = new int[maxvar + 1];
    int *q = new int[maxvar + 1];
    for (int i = 1; i <= maxvar; i++) size[i] = used[i] = 0;
    for (int i = 1; i <= maxvar; i++) {
        int l = 1, r = 0;
        for (int j = 1; j <= maxvar; j++) used[j] = 0;
        q[++r] = i, used[i] = true;
        while (l <= r) {
            int u = q[l++];
            auto& gt = gate[cell[u].gate];
            if(gt.ins == 0) continue;
            assert(gt.ins == 2);
            int a = gt[0], b = gt[1];
            a = abs(a), b = abs(b);
            if(!used[a]) q[++r] = a, used[a] = true;
            if(!used[b]) q[++r] = b, used[b] = true;
        }   
        size[i] = r;
    }
    delete q;
}

// x, y -> z, z, a -> b; z = 0

int Sweep::merge(int x, int y) {
    assert(x > 0 && y > 0 && active[x] && active[y]);
    if (topo_index[x] > topo_index[y]) std::swap(x, y);
    ++stamp;
    for (int i = 0; i < inv_C[x].size(); i++)
        var_stamp[inv_C[x][i]] = stamp;
    for (int i = 0; i < inv_C[y].size(); i++) {
        int c = inv_C[y][i];
        int v = abs(gate[c].out);
        if (!active[v]) continue;
        if (gate[c].ins == 1) {
            int l1 = gate[c][0];
            assert(abs(l1) == y);
            gate[c][0] = gate[c][0] > 0 ? x : -x;
            fanout[x]++;
            inv_C[x].push(c);
            continue;
        }
        int l1 = gate[c][0], l2 = gate[c][1];
        int v1 = abs(l1), v2 = abs(l2);
        assert(active[v1] && active[v2]);
        if (v1 == y) gate[c][0] = gate[c][0] > 0 ? x : -x;
        if (v2 == y) gate[c][1] = gate[c][1] > 0 ? x : -x;
        if (var_stamp[c] == stamp) {
            assert(abs(gate[c][0]) == abs(gate[c][1]));
            if (gate[c][0] == gate[c][1]) {
                gate[c].ins = 1;
                //need simplify
            }
            else {
                // printf("\nc wrong\n");
                // exit(-1);
                // fix_var(gate[c].out);
                fanout[x]++;
                inv_C[x].push(c);
            }
        }
        else {
            fanout[x]++;
            inv_C[x].push(c);
        }
    }
    // return 0;
    int l = 1, r = 0;
    q[++r] = y;
    while (l <= r) {
        int u = q[l++];
        inv_C[u].clear(), fanout[u] = 0, active[u] = 0;
        int c = cell[u].gate;
        for (int i = 0; i < gate[c].ins; i++) {
            int v = abs(gate[c][i]);
            if (!active[v]) continue;
            if (!(--fanout[v])) q[++r] = v;
        }
    }
    return r;
}
// y, t -> z
// x, t -> z
void Sweep::cal_topo_index() {
    for (int i = 1; i <= maxvar; i++) topo_counter[i] = 0;
    int cnt = 0;
    std::queue<int> q;
    for(int i=0; i<inputs; i++) {
        q.push(input[i]);
        topo_index[input[i]] = ++cnt;
    }
    while(!q.empty()) {
        int u = q.front();
        q.pop();
        for(int i=0; i<inv_V[u].size(); i++) {
            int v = inv_V[u][i];
            topo_counter[v]++;
            if(topo_counter[v] == 2) {
                q.push(v);
                topo_index[v] = ++cnt;
            }
        }
    }
}

bool Sweep::match_xor(int x) { // check clause x, -y, -z <==> xor
    int y = gate[x][0], z = gate[x][1];
    if (y >= 0 || z >= 0) return false;
    int gy = cell[-y].gate, gz = cell[-z].gate;
    if (!gy || !gz || gate[gy].type || gate[gy].type) return false;
    int u1 = gate[gy][0], v1 = gate[gy][1];
    int u2 = gate[gz][0], v2 = gate[gz][1];
    if (u1 == -u2 && v1 == -v2) goto doit;
    if (u1 == -v2 && v1 == -u2) goto doit;
    return false;
doit:
    if (cell[-y].fanouts <= 1) cell[-y].gate = 0, gate[gy].type = UNUSE;
    if (cell[-z].fanouts <= 1) cell[-z].gate = 0, gate[gz].type = UNUSE;
    gate[x][0] = gate[gy][0], gate[x][1] = gate[gy][1];
    gate[x].type = Xor;

    return true;
}

void Sweep::recalculate_outs() {
    for (int i = 1; i <= maxvar; i++) 
        cell[i].fanouts = 0;
    for (int i = 1; i < gate.size(); i++) {
        if (gate[i].type == UNUSE) continue;
        for (int j = 0; j < gate[i].ins; j++) {
            cell[abs(gate[i][j])].fanouts++;
        }
    }
}

inline int Sweep::find_par(int x) {
    return x == par[x] ? x : (par[x] = find_par(par[x]));
}

void Sweep::cal_xor() {
    for (int i = 1; i <= maxvar; i++) par[i] = i;
    for (int i = 1; i < gate.size(); i++) {
        if (gate[i].type != Xor) continue;
        for (int j = 0; j < gate[i].ins; j++) {
            int par0 = find_par(abs(gate[i].out));
            par[find_par(abs(gate[i][0]))] = par0;
            par[find_par(abs(gate[i][1]))] = par0;
        }
    }
    int *sign = new int[maxvar + 1];
    for (int i = 1; i <= maxvar; i++) sign[i] = 0;
    for (int i = 1; i <= maxvar; i++) {
        if (i == find_par(i)) continue;
        sign[par[i]] = 1;
    }
    for (int i = 1; i <= maxvar; i++)
        if (sign[i]) sign[i] = ++xors;
    for (int i = 1; i <= maxvar; i++) {
        xors_id[i] = sign[par[i]];
    }
    // printf("c total xors block: %d\n", xors);
}

void Sweep::show_struct() {
    printf("===========================struct==========================\n");
    for (int i = 1; i <= maxvar; i++) {
        int c = cell[i].gate;
        if (!c || gate[c].type == UNUSE) continue;
       
        printf("%d", gate[c].out);

        printf(" = %s ( ", gate_type[gate[c].type].c_str());

        for (int j = 0; j < gate[c].ins; j++) {
            printf("%d ", gate[c][j]);
        }
        puts(")");
    }
    printf("===========================================================\n");
}

void Sweep::identify() {

    recalculate_outs();
    for (int i = 1; i <= maxvar; i++) {
        int c = cell[i].gate;
        if (!c || gate[c].type) continue;
        if (match_xor(c)) continue;
    }
    recalculate_outs();

    
    // cell_cp = new type_cell[maxvar + 1];
    // gate_cp.push();
    // for (int i = 1; i < gate.size(); i++)
    //     gate_cp.push(gate[i]);
    // for (int i = 1; i <= maxvar; i++) 
    //     cell_cp[i] = cell[i];
    

    // show_struct();
    // cal_xor();
 
    inv_C = new vec<int>[maxvar + 1];
    for (int i = 1; i < gate.size(); i++) {
        if (gate[i].type == UNUSE) continue;
        active[abs(gate[i].out)] = 1;
        for (int j = 0; j < gate[i].ins; j++) {
            inv_C[abs(gate[i][j])].push(i), active[abs(gate[i][j])] = 1;
        }
    }

    inv_V = new vec<int>[maxvar + 1];
    for(int i=1; i<=maxvar; i++) {
        int c = cell[i].gate;
        if(gate[c].type == UNUSE) continue;
        for(int j=0; j<gate[c].ins; j++) {
            inv_V[abs(gate[c][j])].push(i);
        }
    }

    for (int i = 1; i <= maxvar; i++) fanout[i] = inv_C[i].size();
}