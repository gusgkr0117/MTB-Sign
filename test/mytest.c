#include <pari/pari.h>
#include <time.h>
#include "ideal.h"
#include "toolbox.h"
#include "klpt.h"
#include "printing.h"

#include "curve.h"
#include "mont.h"
#include "idiso.h"
#include "constants.h"
#include "precomputed.h"

// Evaluating an isogeny
#include "isogenies.h"


struct quaternion_setup_t {
    GEN p; // the prime
    uintbig p_uint; // the prime
    GEN B; // the quaternion algebra
    GEN qf; // the quaternion algebra
    GEN O0; // the cannonical maximal order
    GEN one;
    GEN i;
    GEN j;
    GEN ji;
    GEN torsion_fm; // factorisation matrix of the available torsion

    GEN O0_b1;
    GEN O0_b2;
    GEN O0_b3;
    GEN O0_b4;
    GEN O0_to_standard;
    GEN standard_to_O0;

    proj E0;
};


struct quaternion_setup_t precomp_test;

void random_point(proj *P, proj const *A, long ell, long e, bool twist) {
    uintbig cofactor;
    if ((!twist && curve_order_is_p_plus_one) || (twist && !curve_order_is_p_plus_one)) 
        uintbig_add3(&cofactor, &p, &uintbig_1);
    else { uintbig_sub3(&cofactor, &p, &uintbig_1); }


    uintbig ell_big;
    uintbig_set(&ell_big, ell);
    for (int i = 0; i < e; ++i) {
        uintbig_div3_64(&cofactor, &cofactor, ell); 
    }
    proj Z;




    while (1) {
        fp2_random(&P->x); P->z = fp2_1;
        if (twist == is_on_curve(P, A)) continue;
        xMUL(P, A, P, &cofactor);
        Z = *P;
        for (int i = 0; i < e-1; ++i) {
            xMUL(&Z, A, &Z, &ell_big);
        }
        if (!fp2_iszero(&Z.z)) { 
            // a final test
            xMUL(&Z, A, &Z, &ell_big);
            assert(fp2_iszero(&Z.z));
            return;
        }
    }
}

void compute_action(GEN *M, void (*endo)(point*, const point*), const point *P1, const point *P2, const proj *E, long ell, long e) {
    point Q;
    GEN a,b,c,d;
    assert(ted_is_on_curve(P1,E));
    endo(&Q,P1);
    // proj E1;
    // printf("E0 = %s\n", proj_code(E));
    // fp2_frob2(&E1.x, &E->x);
    // fp2_frob2(&E1.z, &E->z);
    assert(ted_is_on_curve(&Q,E));

    assert(ted_bidim_log(&a, &c, E, &Q, P1, P2, ell, e));

    endo(&Q,P2);
    assert(ted_bidim_log(&b, &d, E, &Q, P1, P2, ell, e));

    *M = mkmat2(mkcol2(a,c),mkcol2(b,d));
}


void random_basis(proj *P1, proj *P2, point *P1_ted, point *P2_ted, proj const *A, long ell, long e, bool twist) {
    point P1_mul, P2_mul, tmp;
    uintbig ell_big;
    uintbig_set(&ell_big, ell);
    fp2 weil;




    proj E;
    mont_to_ted(&E, A, twist);

    random_point(P1, A, ell, e, twist);

    mont_to_ted_point(P1_ted, A, P1);
    P1_mul = *P1_ted;
    for (int i = 0; i < e-1; ++i) {
        ted_mul(&P1_mul, &P1_mul, &E, &ell_big);
    }

    assert(ted_is_on_curve(&P1_mul,&E));
    assert(!ted_iszero(&P1_mul));
    ted_mul(&tmp, &P1_mul, &E, &ell_big);
    assert(ted_iszero(&tmp));

    do {
        random_point(P2, A, ell, e, twist);
        mont_to_ted_point(P2_ted, A, P2);
        P2_mul = *P2_ted;
        for (int i = 0; i < e-1; ++i) {
            ted_mul(&P2_mul, &P2_mul, &E, &ell_big);
        }
        assert(ted_is_on_curve(&P2_mul,&E));
        assert(!ted_iszero(&P2_mul));
        ted_mul(&tmp, &P2_mul, &E, &ell_big);
        assert(ted_iszero(&tmp));

        ted_weil(&weil, &E, &P1_mul, &P2_mul, &ell_big);
        fp2_sub2(&weil, &fp2_1);
    } while (fp2_iszero(&weil));

}



int main() {
    pari_init(80000000, 1<<18);
    init_precomputations();

    setrand(stoi(time(NULL)));
    srand48(time(NULL));

#pragma region[test1]
    // long ell, e;
    // GEN I, v = mkcol2s(0, 0), w, gelle;

    // ell = 84719; e = 12;
    // gelle = powuu(ell, e);

    // gel(v, 1) = randomi(gelle);
    // gel(v, 2) = randomi(gelle);

    // pari_printf("Random kernel basis : %Ps\n", v);

    // I = kernel_to_ideal_O0_ell(v, ell);
    // w = ideal_to_kernel_O0_ell(I, ell);

    // pari_printf("Ideal basis : %Ps\n", lideal_basis(I));
    // pari_printf("Kernel basis : %Ps\n", w);

    // pari_printf("Det(v|w) mod l^e : %Ps\n", gmod(QM_det(mkmat2(v, w)), gelle));
#pragma endregion[test1]

    // Parameter Settings
    long ell_1 = 84719, ell_2 = 3;
    long e_1 = 12, e_2 = 24;
    GEN N_1, N_2;
    proj E, E_twist;
    GEN action_twist_test_2, action_twist_test_3, action_twist_test_4;
    GEN m1 = mkmat2(mkcol2s(1,0),mkcol2s(0,1));
    long q = 3, c = -1;

    init_curve(&precomp_test.E0);
    mont_to_ted(&E, &precomp_test.E0, false);
    mont_to_ted(&E_twist, &precomp_test.E0, true);

    N_1 = powuu(ell_1, e_1);
    N_2 = powuu(ell_2, e_2);

    printf("#Parameter Setting succeeded\n");
    // pick random basis for N_1
    proj basis_twist_N_1[3]; 
    point basis_twist_ted_N_1[3];

    random_basis(&basis_twist_N_1[0], &basis_twist_N_1[1], &basis_twist_ted_N_1[0], &basis_twist_ted_N_1[1], &global_setup.E0, ell_1, e_1, true);

    printf("#Generating Random Basis succeeded\n");

    // Compute action
    GEN m_frob, m_dist;

    compute_action(&m_dist, ted0_dist, &basis_twist_ted_N_1[0], &basis_twist_ted_N_1[1], &E_twist, ell_1, e_1);
    compute_action(&m_frob, ted0_frob_twist, &basis_twist_ted_N_1[0], &basis_twist_ted_N_1[1], &E_twist, ell_1, e_1);

    printf("#Computing Action succeeded\n");

    // Compute actions of each basis element(2-by-2 matrices)
    action_twist_test_2 = gmod(gneg(m_dist),N_1); //(1+i)/2
    action_twist_test_3 = gmod(gmul(action_twist_test_2,m_frob),N_1); //(j-ji)/2
    action_twist_test_4 = gmod(gadd(action_twist_test_2,gmul(action_twist_test_3,stoi(c))),N_1); //(1+i+cj-cji)/2
    action_twist_test_4 = gmod(gmul(action_twist_test_4,gen_2),N_1); // (1+i+cj-cji)
    action_twist_test_4 = gmod(gsub(gsub(action_twist_test_4,gmul(m_frob,stoi(c))),m1),N_1); // (i-cji)

    // pari_printf("%Ps\n", action_twist_test_2);
    // pari_printf("%Ps\n", action_twist_test_3);
    // pari_printf("%Ps\n", action_twist_test_4);

    pari_printf("%Ps\n", global_setup.action_twist_2[3]);
    pari_printf("%Ps\n", global_setup.action_twist_3[3]);
    pari_printf("%Ps\n", global_setup.action_twist_4[3]);

    // Generating random N2-degree isogeny and Computing the codomain curve E_1
    GEN v = mkcol2s(0, 0), random_I;
    odd_isogeny random_phi;
    proj codomain;
    proj kernel_point, Q;

    gel(v, 1) = randomi(N_2);
    gel(v, 2) = randomi(N_2);

    kernel_point = vec_to_E0(v, false);

    random_I = kernel_to_ideal_O0_ell(v, ell_2);

    pari_printf("%Ps\n", lideal_basis(random_I));

    init_curve(&codomain);
    xISOG(&codomain, &Q, &kernel_point, ell_2);

    printf("codomain coeff : \n");
    print_fp2(&codomain.x);
    print_fp2(&codomain.z);

    // Pick random basis on the codomain curve
    proj basis[2];
    point basis_ted[2];
    random_basis(&basis[0], &basis[1], &basis_ted[0], &basis_ted[1], &codomain, ell_1, e_1, true);

    // Computing actions of each basis element of Hom(E_0, E_1)
    GEN action_twist_hom_2, action_twist_hom_3, action_twist_hom_4;


    pari_close();
    return 0;
}