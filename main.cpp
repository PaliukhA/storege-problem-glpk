#include <stdio.h>
#include <stdlib.h>
#include "glpk.h"
#include <vector>
#include <iostream>
#include <cstdio>
#include <iomanip>

using namespace std;
double eps = 1e-4;

int main() {
    cout << setprecision(6);
    cout << fixed;
    int cnt_storage, cnt_client;
//    std::freopen("../capacitated-warehouse-location/01", "r", stdin);
    std::cin >> cnt_storage >> cnt_client;
    std::vector<double> capacity(cnt_storage);
    std::vector<double> open_cost(cnt_storage);
    std::vector<double> demend(cnt_client);
    glp_term_out(GLP_OFF);
    glp_prob *mip;
    mip = glp_create_prob();
    glp_set_obj_dir(mip, GLP_MIN);
    glp_add_rows(mip, cnt_client + cnt_storage);
    for (int i = 0; i < cnt_client; ++i) {
        glp_set_row_bnds(mip, i + 1, GLP_FX, 1, 1);
    }
    for (int i = 0; i < cnt_storage; ++i) {
        glp_set_row_bnds(mip, i + 1 + cnt_client, GLP_LO, 0, 0);
    }
    glp_add_cols(mip, cnt_storage + cnt_storage * cnt_client);
    vector<int> idx_i(1, 0);
    vector<int> idx_j(1, 0);
    vector<double> vals(1, 0);
    for (int i = 0; i < cnt_storage; ++i) {
        for (int j = 0; j < cnt_client; ++j) {
            idx_i.push_back(j + 1);
            idx_j.push_back(cnt_storage + i * cnt_client + j + 1);
            vals.push_back(1);
        }
    }
    for (int i = 0; i < cnt_storage; ++i) {
        std::cin >> capacity[i];
        std::cin >> open_cost[i];
        glp_set_obj_coef(mip, i + 1, open_cost[i]);
        glp_set_col_kind(mip, i + 1, GLP_IV);
        glp_set_col_bnds(mip, i + 1, GLP_DB, 0, 1);

    }
    for (int i = 0; i < cnt_storage; ++i) {
        idx_i.push_back(i + 1 + cnt_client);
        idx_j.push_back(i + 1);
        vals.push_back(capacity[i]);
    }
    for (int i = 0; i < cnt_client; ++i) {
        std::cin >> demend[i];
    }
    std::vector<std::vector<double> > use_cost(cnt_storage);
    for (int i = 0; i < cnt_storage; ++i) {
        for (int j = 0; j < cnt_client; ++j) {
            idx_i.push_back(i + 1 + cnt_client);
            idx_j.push_back(cnt_storage + i * cnt_client + j + 1);
            vals.push_back(-demend[j]);
        }
    }
    for (int i = 0; i < cnt_storage; ++i) {
        for (int j = 0; j < cnt_client; ++j) {
            double val;
            std::cin >> val;
            use_cost[i].push_back(val);
            int index = cnt_storage + i * cnt_client + j + 1;
            glp_set_obj_coef(mip, index, val);
            glp_set_col_bnds(mip, index, GLP_LO, 0, 0);
        }
    }
    glp_load_matrix(mip, static_cast<int>(vals.size()) - 1, idx_i.data(), idx_j.data(), vals.data());
    glp_smcp params_lp;
    glp_init_smcp(&params_lp);
    params_lp.msg_lev = GLP_MSG_OFF;
    glp_simplex(mip, &params_lp);
    glp_iocp param_mip;
    glp_init_iocp(&param_mip);
    param_mip.msg_lev = GLP_MSG_OFF;
    param_mip.br_tech = GLP_BR_PCH;
    param_mip.gmi_cuts = GLP_ON;
    param_mip.mir_cuts = GLP_ON;
    param_mip.presolve = GLP_ON;
    param_mip.binarize = GLP_ON;

    glp_intopt(mip, &param_mip);
    vector<int> take(cnt_storage, 0);
    int cnt_ans = 0;
    for (int i = 0; i < cnt_storage; ++i) {
        if (glp_mip_col_val(mip, i + 1) > eps) {
            take[i] = 1;
            cnt_ans++;
        }
    }
    cout << cnt_ans << '\n';
    for (int i = 0; i < cnt_storage; ++i) {
        if(take[i] == 1) {
            cout << i+1 << ' ';
        }
    }
    cout << '\n';
    for (int i = 0; i < cnt_storage; ++i) {
        if (take[i] == 1) {
            for (int j = 0; j < cnt_client; ++j) {
                cout << glp_mip_col_val(mip, cnt_storage + i * cnt_client + j + 1) << ' ';
            }
            cout << '\n';
        }
    }
    cout << '\n';
    glp_delete_prob(mip);
    return 0;
}

