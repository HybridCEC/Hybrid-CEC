#include "sweep.hpp"

extern "C" {
   #include "src/internal.h"
}

bool Sweep::check_var_struct(int x, int y) {
    if (size[x] != size[y]) return false;
    int *use1 = new int[maxvar + 1];
    int *use2 = new int[maxvar + 1];
    for (int i = 1; i <= maxvar; i++) use1[i] = use2[i] = 0;
    std::queue<int> q1, q2;

    q1.push(abs(x)), use1[abs(x)] = true;
    q2.push(abs(y)), use2[abs(y)] = true;
    
    while(!q1.empty() && !q2.empty()) {
        int u = q1.front();
        q1.pop();
        int v = q2.front();
        q2.pop();
        auto& gt1 = gate[cell[u].gate];
        auto& gt2 = gate[cell[v].gate];
        if (gt1.ins != gt2.ins || gt1.type != gt2.type) return false;
        if(gt1.ins == 0) continue;
        int a = gt1[0], b = gt1[1];
        int c = gt2[0], d = gt2[1];
        int va = abs(a), vb = abs(b), vc = abs(c), vd = abs(d);
        if (size[va] != size[vc]) std::swap(c, d), std::swap(vc, vd);
        if (size[va] != size[vc] || size[vb] != size[vd]) return false;
        if (size[va] == size[vb]) {
            if (sign(a) != sign(c)) std::swap(c, d), std::swap(vc, vd);
            if (use1[va] != use2[vc]) std::swap(c, d), std::swap(vc, vd);
        }
        if (sign(a) != sign(c) || sign(b) != sign(d)) return false;
        if (use1[va] != use2[vc] || use1[vb] != use2[vd]) return false;
        if(!use1[va]) q1.push(va), q2.push(vc), use1[va] = use2[vc] = true;
        if(!use1[vb]) q1.push(vb), q2.push(vd), use1[vb] = use2[vd] = true;
    }
    if (!q1.empty() || !q2.empty()) return false;
    return true;
}

bool Sweep::seq_check_var_equal(int x, int y) {

    ksolver = kissat_init();
    memset(used, 0, sizeof(int) * (maxvar + 1));
    int l = 1, r = 0, connect = 0;
    q[++r] = abs(x), used[abs(x)] = true;
    q[++r] = abs(y), used[abs(y)] = true;
    ++stamp;
    int s = 0, is = 0;
    epcec_in.clear();
    while(l <= r) {
        int u = q[l++];
        var_stamp[u] = stamp;
        if (is_input[u]) epcec_in.push(u), ++is;
        if (xors_id[u]) xors_bot[xors_id[u]]++;
        s++;
        auto& gt = gate[cell[u].gate];
        if(gt.ins == 0) continue;
        if (gt.ins == 1) {
            int a = gt[0];
            kissat_add(ksolver, u), kissat_add(ksolver, -a); kissat_add(ksolver, 0);
            kissat_add(ksolver, -u), kissat_add(ksolver, a); kissat_add(ksolver, 0);
            a = abs(a);
            if (!used[a]) q[++r] = a, used[a] = true;
            continue;
        }
        int a = gt[0], b = gt[1];
        if (gt.type == And) {
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
        assert(active[a] && active[b]);
        if(!used[a]) q[++r] = a, used[a] = true;
        if(!used[b]) q[++r] = b, used[b] = true;
    }
    for (int i = 1; i <= r; i++) par[q[i]] = q[i], xors_bot[q[i]] = 0;
    for (int i = 1; i <= r; i++) {
        int u = q[i], c = cell[u].gate;
        if (is_input[u]) continue;
        if (gate[c].type == Xor) {
            int a = abs(gate[c][0]), b = abs(gate[c][1]);
            int c1 = cell[a].gate, c2 = cell[b].gate;
            if (!c1 || !c2 || gate[c1].type != Xor || gate[c2].type != Xor) continue;
            par[find_par(a)] = find_par(u);
            par[find_par(b)] = find_par(u);
        }
        else if (gate[c].type == And) {
            if (is_input[abs(gate[c][0])] && is_input[abs(gate[c][1])]) ++connect;
        }
        // if (find_par(x) == find_par(y)) printf("c %d %d %d %d cs\n", i, u, a, b);
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
    double score = mx;
    for (int i = 3; i < mx; i++) {
        bucket[i + 1] += bucket[i] >> 2;
        bucket[i] %= 2;
    }
    while (bucket[mx] >= 4) score++, bucket[mx] /= 4;
    for (int i = 1; i <= mx; i++) bucket[i] = 0;
    score /= 1.0 * is;
    double con = 2.0 * connect / (is * (is - 1));
    if (is > 20) {
        max_score = std::max(max_score, score);
        max_connect = std::max(max_connect, con);
    }
    // printf("c score is %.2lf, connect is %.2lf\n", score, con);
    // printf("c inputs number: %d\n", is);
    // if ((is <= 20) || (is <= 32 && score >= 0.3 && con >= 0.3)) {
    if (is <= 32 && score >= 0.3) {
        this_epcec = 1;
        // printf("c USE epcec\n");
        delete ksolver;
        epcec_out.clear();
        epcec_out.push(x);
        epcec_out.push(y);
        return epcec(1); 
    }
    this_epcec = 2;
    // for (int i = 0; i < verified_pairs.size(); i++) {
    //     int x = verified_pairs[i].first, y = verified_pairs[i].second;
    //     kissat_add(ksolver, -x); kissat_add(ksolver, y); kissat_add(ksolver, 0);
    //     kissat_add(ksolver, x); kissat_add(ksolver, -y); kissat_add(ksolver, 0);
    // }

    kissat_add(ksolver, -x), kissat_add(ksolver, -y), kissat_add(ksolver, 0);
    kissat_add(ksolver, x), kissat_add(ksolver, y), kissat_add(ksolver, 0);
    int result = kissat_solve(ksolver);  
    // printf("c var-equal-checking (%d, %d): %d\n", x, y, result);
    delete ksolver;
    return result == 20;
}


bool Sweep::seq_miter() {
    ksolver = kissat_init();
    memset(used, 0, sizeof(int) * (maxvar + 1));
    std::queue<int> q;
    int o = output >> 1;
    q.push(o), used[o] = true;
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
        assert(active[a] && active[b]);
        if(!used[a]) q.push(a), used[a] = true;
        if(!used[b]) q.push(b), used[b] = true;
    }    
    // for (int i = 0; i < verified_pairs.size(); i++) {
    //     int x = verified_pairs[i].first, y = verified_pairs[i].second;
    //     kissat_add(ksolver, -x); kissat_add(ksolver, y); kissat_add(ksolver, 0);
    //     kissat_add(ksolver, x); kissat_add(ksolver, -y); kissat_add(ksolver, 0);
    // }
    if(output & 1) {
        kissat_add(ksolver, -o); kissat_add(ksolver, 0);
    } else {
        kissat_add(ksolver, o); kissat_add(ksolver, 0);
    }
    int result = kissat_solve(ksolver);  

    if(result == 10) {
        int input_vector[inputs];
        for(int i=0; i<inputs; i++) {
            int kval = kissat_value(ksolver, input[i]);
            if(kval > 0) {
                input_vector[i] = 1;
            } else if(kval < 0) {
                input_vector[i] = 0;
            } else {
                input_vector[i] = rand() % 2;
            }
        }
        is_sat = true;
        print_result_and_exit("c find answer by miter!", input_vector);
    }
    return result == 20;
}