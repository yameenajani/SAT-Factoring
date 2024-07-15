/***************************************************************************************[Solver.cc]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include <cstdio>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <stdio.h>
#include <stdarg.h>
#include <obstack.h>
#include <stdint.h>
#include <gmp.h>
#include <cstring>
#include <fplll.h>
//#include "flint/fmpz_polyxx.h"
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz_poly_factor.h"


#include <math.h>
#include <array>
#include <vector>
#include <time.h>

#include "mtl/Sort.h"
#include "core/Solver.h"
#include "utils/System.h"
#include "utils/Options.h"

using namespace Minisat;
using namespace std;
using namespace fplll;
//using namespace flint;

//=================================================================================================
// Options:

void printClause(vec<Lit>& clause) {
    for(int i = 0; i < clause.size(); i++) {
        printf("%s%d ", sign(clause[i]) ? "-" : "", var(clause[i])+1);
    }
    printf("0\n");
    fflush(stdout);
}

static const char* _cat = "CORE";

#if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB
static DoubleOption  opt_step_size         (_cat, "step-size",   "Initial step size",                             0.40,     DoubleRange(0, false, 1, false));
static DoubleOption  opt_step_size_dec     (_cat, "step-size-dec","Step size decrement",                          0.000001, DoubleRange(0, false, 1, false));
static DoubleOption  opt_min_step_size     (_cat, "min-step-size","Minimal step size",                            0.06,     DoubleRange(0, false, 1, false));
#endif
#if BRANCHING_HEURISTIC == VSIDS
static DoubleOption  opt_var_decay         (_cat, "var-decay",   "The variable activity decay factor",            0.95,     DoubleRange(0, false, 1, false));
#endif
#if ! LBD_BASED_CLAUSE_DELETION
static DoubleOption  opt_clause_decay      (_cat, "cla-decay",   "The clause activity decay factor",              0.999,    DoubleRange(0, false, 1, false));
#endif
static DoubleOption  opt_random_var_freq   (_cat, "rnd-freq",    "The frequency with which the decision heuristic tries to choose a random variable", 0, DoubleRange(0, true, 1, true));
static DoubleOption  opt_random_seed       (_cat, "rnd-seed",    "Used by the random variable selection",         91648253, DoubleRange(0, false, HUGE_VAL, false));
static IntOption     opt_ccmin_mode        (_cat, "ccmin-mode",  "Controls conflict clause minimization (0=none, 1=basic, 2=deep)", 2, IntRange(0, 2));
static IntOption     opt_phase_saving      (_cat, "phase-saving", "Controls the level of phase saving (0=none, 1=limited, 2=full)", 2, IntRange(0, 2));
static BoolOption    opt_rnd_init_act      (_cat, "rnd-init",    "Randomize the initial activity", false);
static BoolOption    opt_luby_restart      (_cat, "luby",        "Use the Luby restart sequence", true);
static IntOption     opt_restart_first     (_cat, "rfirst",      "The base restart interval", 100, IntRange(1, INT32_MAX));
static DoubleOption  opt_restart_inc       (_cat, "rinc",        "Restart interval increase factor", 2, DoubleRange(1, false, HUGE_VAL, false));
static DoubleOption  opt_garbage_frac      (_cat, "gc-frac",     "The fraction of wasted memory allowed before a garbage collection is triggered",  0.20, DoubleRange(0, false, HUGE_VAL, false));
static BoolOption    opt_only_sat          (_cat, "only-sat",    "Run only SAT without callback", false);
static BoolOption    opt_hi                (_cat, "hi",          "Call coppersmith with high bits", false);
static BoolOption    opt_lo                (_cat, "lo",          "Call coppersmith with low bits", true);
static BoolOption    opt_cs_both_primes    (_cat, "cs-both-primes", "Call coppersmith with low bits", false);
static BoolOption    opt_prog_hc           (_cat, "prog-hc",     "Use programmatic Heninger/Shacham constraints for branching", false);
static IntOption     opt_cb_wait           (_cat, "cb-wait",     "How many conflicts to wait before calling Coppersmith", 1, IntRange(1, 10000));
// static StringOption  opt_cb_num            (_cat, "cb-num",      "The co-prime to be factored");
// static IntOption     opt_p_lsb_var         (_cat, "p-lsb-var",   "The variable correseponding to the LSB of p", 1, IntRange(1, INT32_MAX));
// static IntOption     opt_p_msb_var         (_cat, "p-msb-var",   "The variable correseponding to the MSB of p", 1, IntRange(1, INT32_MAX));
// static IntOption     opt_q_lsb_var         (_cat, "q-lsb-var",   "The variable correseponding to the LSB of q", 1, IntRange(1, INT32_MAX));
// static IntOption     opt_q_msb_var         (_cat, "q-msb-var",   "The variable correseponding to the MSB of q", 1, IntRange(1, INT32_MAX));
#if BRANCHING_HEURISTIC == CHB
static DoubleOption  opt_reward_multiplier (_cat, "reward-multiplier", "Reward multiplier", 0.9, DoubleRange(0, true, 1, true));
#endif

mpz_t mpz_NUM;
mpz_t mpz_X, pq_temp, p_temp, q_temp, pow2, pow4, inv_factor;
//mpz_t cs_prime;
int prime_len;
int p_msb_count;
vector<char> bits_msb_p1;
vector<char> bits_lsb_p1;
vector<char> bits_msb_p2;
vector<char> bits_lsb_p2;

//=================================================================================================
// Constructor/Destructor:


Solver::Solver() :

    // Parameters (user settable):
    //
    verbosity        (0)
#if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB
  , step_size        (opt_step_size)
  , step_size_dec    (opt_step_size_dec)
  , min_step_size    (opt_min_step_size)
#endif
#if BRANCHING_HEURISTIC == VSIDS
  , var_decay        (opt_var_decay)
#endif
#if ! LBD_BASED_CLAUSE_DELETION
  , clause_decay     (opt_clause_decay)
#endif
  , random_var_freq  (opt_random_var_freq)
  , random_seed      (opt_random_seed)
  , luby_restart     (opt_luby_restart)
  , ccmin_mode       (opt_ccmin_mode)
  , phase_saving     (opt_phase_saving)
  , rnd_pol          (false)
  , rnd_init_act     (opt_rnd_init_act)
  , garbage_frac     (opt_garbage_frac)
  , restart_first    (opt_restart_first)
  , restart_inc      (opt_restart_inc)
  , clbk_st          (opt_only_sat)
  , clbk_hi_bits     (opt_hi)
  , clbk_lo_bits     (opt_lo)
  , cs_both_primes   (opt_cs_both_primes)
  , prog_hc          (opt_prog_hc)
  , cb_wait          (opt_cb_wait)
//   , cb_num           (opt_cb_num)
//   , p_lsb_var        (opt_p_lsb_var)
//   , p_msb_var        (opt_p_msb_var)
//   , q_lsb_var        (opt_q_lsb_var)
//   , q_msb_var        (opt_q_msb_var)

    // Parameters (the rest):
    //
  , learntsize_factor((double)1/(double)3), learntsize_inc(1.1)

    // Parameters (experimental):
    //
  , learntsize_adjust_start_confl (100)
  , learntsize_adjust_inc         (1.5)

    // Statistics: (formerly in 'SolverStats')
    //
  , solves(0), starts(0), decisions(0), rnd_decisions(0), propagations(0), conflicts(0)
  , dec_vars(0), clauses_literals(0), learnts_literals(0), max_literals(0), tot_literals(0)

  , lbd_calls(0)
#if BRANCHING_HEURISTIC == CHB
  , action(0)
  , reward_multiplier(opt_reward_multiplier)
#endif

  , ok                 (true)
#if ! LBD_BASED_CLAUSE_DELETION
  , cla_inc            (1)
#endif
#if BRANCHING_HEURISTIC == VSIDS
  , var_inc            (1)
#endif
  , watches            (WatcherDeleted(ca))
  , qhead              (0)
  , simpDB_assigns     (-1)
  , simpDB_props       (0)
  , order_heap         (VarOrderLt(activity))
  , progress_estimate  (0)
  , remove_satisfied   (true)

    // Resource constraints:
    //
  , conflict_budget    (-1)
  , propagation_budget (-1)
  , asynch_interrupt   (false)
{}

// Global Variables
bool val;
bool call_coppersmith_msb_p1 = true;
bool call_coppersmith_lsb_p1 = true;
bool call_coppersmith_msb_p2 = true;
bool call_coppersmith_lsb_p2 = true;
long int cs_count = 0;
long int cs_count_msb_p1 = 0;
long int cs_count_lsb_p1 = 0;
long int cs_count_msb_p2 = 0;
long int cs_count_lsb_p2 = 0;
long int hc_count = 0;
long int total_hc = 0;
double cs_time = 0;
double hc_time = 0;


Solver::~Solver()
{
    mpz_clear(mpz_NUM);
    //mpz_clear(cs_prime);
    mpz_clears(mpz_X, pq_temp, p_temp, q_temp, pow2, pow4, inv_factor, NULL);
}

void Solver::setInstanceVariables()
{
    mpz_init_set_str(mpz_NUM, cb_num.c_str(), 10);
    //mpz_init(cs_prime);
    mpz_inits(mpz_X, pq_temp, p_temp, q_temp, pow2, pow4, inv_factor, NULL);
    mpz_root(mpz_X, mpz_NUM, 5);
    mpz_tdiv_q_2exp(mpz_X, mpz_X, 2);
    prime_len = p_msb_var;
    p_msb_count = prime_len - (mpz_sizeinbase(mpz_X, 2) - 1); // Subtract 1 as mpz_sizeinbase(mpz_X, 2) gives floor(lg(X))+1
    mpz_set_ui(inv_factor, 1);
    mpz_mul_2exp(inv_factor, inv_factor, p_msb_count);
    mpz_invert(inv_factor, inv_factor, mpz_NUM);
    bits_msb_p1.reserve(prime_len);
    bits_lsb_p1.reserve(p_msb_count);
    bits_msb_p2.reserve(prime_len);
    bits_lsb_p2.reserve(p_msb_count);
}

void Solver::printStats()
{
    double cpu_time = cpuTime();
    double mem_used = memUsedPeak();
    printf("restarts              : %"PRIu64"\n", starts);
    printf("conflicts             : %-12"PRIu64"   (%.0f /sec)\n", conflicts   , conflicts   /cpu_time);
    printf("decisions             : %-12"PRIu64"   (%4.2f %% random) (%.0f /sec)\n", decisions, (float)rnd_decisions*100 / (float)decisions, decisions   /cpu_time);
    printf("propagations          : %-12"PRIu64"   (%.0f /sec)\n", propagations, propagations/cpu_time);
    printf("conflict literals     : %-12"PRIu64"   (%4.2f %% deleted)\n", tot_literals, (max_literals - tot_literals)*100 / (double)max_literals);
    if(!opt_only_sat && (opt_lo || opt_hi)) {
        printf("CS Clauses            : %-12"PRIu64"   (p: hi %ld lo %ld; q: hi %ld lo %ld)\n", cs_count, cs_count_msb_p1, cs_count_lsb_p1, cs_count_msb_p2, cs_count_lsb_p2);
        printf("Total CS time         : %g s\n", cs_time);
        printf("Average CS time       : %g s\n", cs_time/cs_count);
    }
    if(!opt_only_sat && opt_prog_hc) {
        printf("HC Prog. Branchings   : %-12"PRIu64"   (%ld recomputations)\n", hc_count, total_hc);
        printf("Total HC time         : %g s\n", hc_time);
    }
    if (mem_used != 0) printf("Memory used           : %.2f MB\n", mem_used);
    printf("CPU time              : %g s\n", cpu_time);
}

//=================================================================================================
// Minor methods:


// Creates a new SAT variable in the solver. If 'decision' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
Var Solver::newVar(bool sign, bool dvar)
{
    int v = nVars();
    watches  .init(mkLit(v, false));
    watches  .init(mkLit(v, true ));
    assigns  .push(l_Undef);
    vardata  .push(mkVarData(CRef_Undef, 0));
    //activity .push(0);
    activity .push(rnd_init_act ? drand(random_seed) * 0.00001 : 0);
    seen     .push(0);
    polarity .push(sign);
    decision .push();
    trail    .capacity(v+1);
    lbd_seen.push(0);
    picked.push(0);
    conflicted.push(0);
#if ALMOST_CONFLICT
    almost_conflicted.push(0);
#endif
#if ANTI_EXPLORATION
    canceled.push(0);
#endif
#if BRANCHING_HEURISTIC == CHB
    last_conflict.push(0);
#endif
    total_actual_rewards.push(0);
    total_actual_count.push(0);
    setDecisionVar(v, dvar);
    return v;
}


bool Solver::addClause_(vec<Lit>& ps)
{
    assert(decisionLevel() == 0);
    if (!ok) return false;

    // Check if clause is satisfied and remove false/duplicate literals:
    sort(ps);
    Lit p; int i, j;
    for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
        if (value(ps[i]) == l_True || ps[i] == ~p)
            return true;
        else if (value(ps[i]) != l_False && ps[i] != p)
            ps[j++] = p = ps[i];
    ps.shrink(i - j);

    if (ps.size() == 0)
        return ok = false;
    else if (ps.size() == 1){
        uncheckedEnqueue(ps[0]);
        return ok = (propagate() == CRef_Undef);
    }else{
        CRef cr = ca.alloc(ps, false);
        clauses.push(cr);
        attachClause(cr);
    }

    return true;
}


void Solver::attachClause(CRef cr) {
    const Clause& c = ca[cr];
    assert(c.size() > 1);
    watches[~c[0]].push(Watcher(cr, c[1]));
    watches[~c[1]].push(Watcher(cr, c[0]));
    if (c.learnt()) learnts_literals += c.size();
    else            clauses_literals += c.size(); }


void Solver::detachClause(CRef cr, bool strict) {
    const Clause& c = ca[cr];
    assert(c.size() > 1);
    
    if (strict){
        remove(watches[~c[0]], Watcher(cr, c[1]));
        remove(watches[~c[1]], Watcher(cr, c[0]));
    }else{
        // Lazy detaching: (NOTE! Must clean all watcher lists before garbage collecting this clause)
        watches.smudge(~c[0]);
        watches.smudge(~c[1]);
    }

    if (c.learnt()) learnts_literals -= c.size();
    else            clauses_literals -= c.size(); }


void Solver::removeClause(CRef cr) {
    Clause& c = ca[cr];
    detachClause(cr);
    // Don't leave pointers to free'd memory!
    if (locked(c)) vardata[var(c[0])].reason = CRef_Undef;
    c.mark(1); 
    ca.free(cr);
}


bool Solver::satisfied(const Clause& c) const {
    for (int i = 0; i < c.size(); i++)
        if (value(c[i]) == l_True)
            return true;
    return false; }

int bit_i = 0;
// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
//
void Solver::cancelUntil(int level) {
    if (decisionLevel() > level){
        for (int c = trail.size()-1; c >= trail_lim[level]; c--){
            Var      x  = var(trail[c]);
            if((x >= 0 && x < bit_i) || (x >= q_msb_var - prime_len && x < q_msb_var - prime_len + bit_i)) {
                // Bit in LSB of prime1 or prime2 being unset, so we will recompute pq_temp, p_temp, q_temp
                bit_i = 0;
            }
            uint64_t age = conflicts - picked[x];
            if (age > 0) {
                double reward = ((double) conflicted[x]) / ((double) age);
#if BRANCHING_HEURISTIC == LRB
#if ALMOST_CONFLICT
                double adjusted_reward = ((double) (conflicted[x] + almost_conflicted[x])) / ((double) age);
#else
                double adjusted_reward = reward;
#endif
                double old_activity = activity[x];
                activity[x] = step_size * adjusted_reward + ((1 - step_size) * old_activity);
                if (order_heap.inHeap(x)) {
                    if (activity[x] > old_activity)
                        order_heap.decrease(x);
                    else
                        order_heap.increase(x);
                }
#endif
                total_actual_rewards[x] += reward;
                total_actual_count[x] ++;
            }
#if ANTI_EXPLORATION
            canceled[x] = conflicts;
#endif
            assigns [x] = l_Undef;
            if (phase_saving > 1 || (phase_saving == 1) && c > trail_lim.last())
                polarity[x] = sign(trail[c]);
            insertVarOrder(x); }
        qhead = trail_lim[level];
        trail.shrink(trail.size() - trail_lim[level]);
        trail_lim.shrink(trail_lim.size() - level);
    } }


//=================================================================================================
// Major methods:


Lit Solver::pickBranchLit()
{
    Var next = var_Undef;

    // Random decision:
    if (drand(random_seed) < random_var_freq && !order_heap.empty()){
        next = order_heap[irand(random_seed,order_heap.size())];
        if (value(next) == l_Undef && decision[next])
            rnd_decisions++; }

    // Activity based decision:
    while (next == var_Undef || value(next) != l_Undef || !decision[next])
        if (order_heap.empty()){
            next = var_Undef;
            break;
        } else {
#if ANTI_EXPLORATION
            next = order_heap[0];
            uint64_t age = conflicts - canceled[next];
            while (age > 0 && value(next) == l_Undef) {
                double decay = pow(0.95, age);
                activity[next] *= decay;
                if (order_heap.inHeap(next)) {
                    order_heap.increase(next);
                }
                canceled[next] = conflicts;
                next = order_heap[0];
                age = conflicts - canceled[next];
            }
#endif
            next = order_heap.removeMin();
        }

    return next == var_Undef ? lit_Undef : mkLit(next, rnd_pol ? drand(random_seed) < 0.5 : polarity[next]);
}

// CUSTOM COPPERSMITH IMPLEMENTATION
void nCr(mpz_t rop, int n, int r, mpz_t p_tilda) {
    mpz_t t1, t2;
    mpz_inits(t1, t2, NULL);

    mpz_bin_uiui(t1, n, r);
    mpz_pow_ui(t2, p_tilda, n-r);
    mpz_mul(rop, t1, t2);

    mpz_clears(t1, t2, NULL);
}

bool coppersmith(mpz_t p_tilda_num, int is_lsb) {
    cs_count += 1;
    mpz_t prime, t, t2, p_tilda, mpz_val;
    unsigned int len_knownbits = p_msb_count, h = 2, k;
    k = 2*h;

    mpz_inits(prime, t, t2, p_tilda, mpz_val, NULL);

    if(is_lsb) {
        mpz_mul(p_tilda, inv_factor, p_tilda_num);
        mpz_mod(p_tilda, p_tilda, mpz_NUM); // Reduce constant term mod N
    } else {
        mpz_set(p_tilda, p_tilda_num);
    }

    mpz_t polys[k+1][k+1];
    for(int i = 0; i < k+1; i++) {
        for(int j = 0; j < k+1; j++) {
            mpz_init(polys[i][j]);
        }
    }

    int z;
    mpz_t coeff;
    mpz_init(coeff);

    for(int i = 0; i < h+1; i++) {
        mpz_pow_ui(t, mpz_NUM, h-i);
        for(int j = 0; j < i+1; j++) {
            nCr(t2, i, j, p_tilda);
            mpz_mul(coeff, t, t2);
            mpz_pow_ui(t2, mpz_X, j);
            mpz_mul(coeff, coeff, t2);
            mpz_set(polys[i][j], coeff);
        }
        z = i;
    }
    z++;
    for(int i = 1; i < k-h+1; i++) {
        for(int j = 0; j < h+1; j++) {
            nCr(t, h, j, p_tilda);
            mpz_set(polys[z][j+i], t);
        }
        for(int j = 0; j < h+i+1; j++) {
            mpz_pow_ui(t2, mpz_X, j);
            mpz_mul(coeff, polys[z][j], t2);
            mpz_set(polys[z][j], coeff);
        }
        z++;
    }

    mpz_clear(coeff);

    ZZ_mat<mpz_t> * FPlllMat;
    Z_NR<mpz_t>  zval;

    FPlllMat = new ZZ_mat<mpz_t>(k+1,k+1);

    for(int i = 0; i < k+1; i++) {
      for(int j = 0; j < k+1; j++) {
        zval = polys[i][j];
        (*FPlllMat)[i][j] = zval;
      }
    }

    lll_reduction(*FPlllMat, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER, FT_DEFAULT, 0, LLL_DEFAULT);

    fmpz_poly_t pol;
    fmpz_poly_init(pol);

    mpz_t red_coef, temp_var;
    mpz_inits(red_coef, temp_var, NULL);

    fmpz_t fmpz_var;
    fmpz_init(fmpz_var);

    for(int i=0; i<k+1; i++) {
        (*FPlllMat)[0][i].get_mpz(red_coef);
        mpz_pow_ui(temp_var, mpz_X, i);
        mpz_div(red_coef, red_coef, temp_var);
        fmpz_set_mpz(fmpz_var, red_coef);
        fmpz_poly_set_coeff_fmpz(pol, i, fmpz_var);
    }

    mpz_clears(temp_var, red_coef, NULL);

    fmpz_poly_factor_t factors;
    fmpz_poly_factor_init(factors);

    fmpz_poly_factor_zassenhaus(factors, pol);

    for(int i = 0; i < factors->num; i++)
    {
        if(fmpz_poly_degree(&factors->p[i]) == 1 && fmpz_poly_get_coeff_ui(&factors->p[i], 1) == 1) {
            fmpz_poly_get_coeff_fmpz(fmpz_var, &factors->p[i], 0);
            fmpz_get_mpz(mpz_val, fmpz_var);
            if(mpz_cmp_ui(mpz_val, 0)!=0)
            {
                mpz_mul_si(mpz_val, mpz_val, -1);
                if(is_lsb) {
                    mpz_mul_2exp(prime, mpz_val, len_knownbits);
                    mpz_add(prime, prime, p_tilda_num);
                } else {
                    mpz_add(prime, p_tilda_num, mpz_val);
                }
                if(mpz_divisible_p(mpz_NUM, prime)) {
                    val = true;
                    //mpz_set(cs_prime, prime);
                    printf("============================================================================================"); if(opt_prog_hc) printf("============="); printf("\n");
                    gmp_printf("The result from Coppersmith is - %Zd\n", mpz_val);
                    gmp_printf("One of the factors of N is - %Zd\n", prime);
                    break;
                } else {
                    val = false;
                }
            }
        } else {
            val = false;
        }
    }

    fmpz_clear(fmpz_var);

    fmpz_poly_factor_clear(factors);
    fmpz_poly_clear(pol);

    for(int i = 0; i < k+1; i++) {
      for(int j = 0; j < k+1; j++) {
        mpz_clear(polys[i][j]);
      }
    }

    mpz_clears(prime, t, t2, p_tilda, mpz_val, NULL);

    delete FPlllMat;

    return val;
}

long next_cs_call = 0;
int wait = 0;
// A callback function for programmatic interface. If the callback detects conflicts, then
// refine the clause database by adding clauses to out_learnts. This function is called
// very frequently, if the analysis is expensive then add code to skip the analysis on
// most calls. However, if complete is set to true, do not skip the analysis or else the
// solver will be unsound.
//
// complete: true if and only if the current trail is a complete assignment that satisfies
//           the clause database. Note that not every variable is necessarily assigned since
//           the simplification steps may have removed some variables! If complete is true,
//           the solver will return satisfiable immediately unless this function returns at
//           least one clause.

void Solver::callbackFunction(bool complete, vec<vec<Lit> >& out_learnts, Lit& branch_lit) {
    if(opt_only_sat) return;

    bool start_coppersmith = true;
    if(opt_prog_hc) {
        double hc_start = cpuTime();
        int plsb = p_lsb_var - 1 + bit_i;
        int qlsb = q_lsb_var - 1 + bit_i;
        int ptrack = plsb;
        int qtrack = qlsb;

        int pbit, qbit;

        if(bit_i == 0) { // Compute pq_temp, p_temp, q_temp from scratch
            mpz_set_ui(pq_temp, 0);
            mpz_set_ui(p_temp, 0);
            mpz_set_ui(q_temp, 0);
            mpz_set_ui(pow2, 1);
            mpz_set_ui(pow4, 1);
            total_hc++; // Total number of recomputions
        }

        while(assigns[plsb]!=l_Undef && assigns[qlsb]!=l_Undef && bit_i < prime_len-1)
        {
            if(assigns[qlsb] == l_True) {
                mpz_addmul(pq_temp, p_temp, pow2);
            }
            if(assigns[plsb] == l_True) {
                mpz_addmul(pq_temp, q_temp, pow2);
                mpz_add(p_temp, p_temp, pow2);
            }
            if(assigns[qlsb] == l_True) {
                mpz_add(q_temp, q_temp, pow2);
            }
            if(assigns[qlsb] == l_True && assigns[plsb] == l_True) {
                mpz_add(pq_temp, pq_temp, pow4);
            }

            mpz_mul_ui(pow2, pow2, 2);
            mpz_mul_ui(pow4, pow4, 4);

            plsb++;
            qlsb++;
            bit_i++;
        }

        mpz_t temp;
        mpz_init(temp);
        mpz_sub(temp, mpz_NUM, pq_temp);
        int bitValue = mpz_tstbit(temp, bit_i);
        mpz_clear(temp);

        // Determine if the Heninger/Shacham constraints imply the value of an unassigned variable;
        // If so, then branch on that variable.
        if(assigns[plsb]==l_Undef && assigns[qlsb]==l_True) {
            branch_lit = mkLit(plsb, bitValue == 0);
        } else if(assigns[plsb]==l_Undef && assigns[qlsb]==l_False) {
            branch_lit = mkLit(plsb, bitValue == 1);
        } else if(assigns[qlsb]==l_Undef && assigns[plsb]==l_True) {
            branch_lit = mkLit(qlsb, bitValue == 0);
        } else if(assigns[qlsb]==l_Undef && assigns[plsb]==l_False) {
            branch_lit = mkLit(qlsb, bitValue == 1);
        }

        hc_time += (cpuTime() - hc_start);

        if(branch_lit != lit_Undef) {
            hc_count++;
            start_coppersmith = false;
        }
    }

    if(start_coppersmith && (wait%opt_cb_wait==0) && next_cs_call <= conflicts)
    {
        next_cs_call = conflicts + 1;
        call_coppersmith_msb_p1 = true;
        call_coppersmith_msb_p2 = true;
        call_coppersmith_lsb_p1 = true;
        call_coppersmith_lsb_p2 = true;
        bits_msb_p1.clear();
        bits_msb_p2.clear();
        bits_lsb_p1.clear();
        bits_lsb_p2.clear();
        
        // Checking if msb bits of p1 equivalent to p_msb_count are set
        for(int i=p_msb_var-1; i>p_msb_var-p_msb_count-1; i--)
        {
            if(assigns[i]==l_Undef) {
                call_coppersmith_msb_p1 = false;
                break;
            } else {
                bits_msb_p1.push_back((assigns[i]==l_True ? '1' : '0'));
            }
        }

        // Checking if msb bits of p2 equivalent to p_msb_count are set
        for(int i=q_msb_var-1; i>q_msb_var-p_msb_count-1; i--)
        {
            if(assigns[i]==l_Undef) {
                call_coppersmith_msb_p2 = false;
                break;
            } else {
                bits_msb_p2.push_back((assigns[i]==l_True ? '1' : '0'));
            }
        }

        // Checking if lsb bits of p1 equivalent to p_msb_count are set
        for(int i=p_lsb_var-1; i<p_lsb_var+p_msb_count-1; i++)
        {
            if(assigns[i]==l_Undef) {
                call_coppersmith_lsb_p1 = false;
                break;
            } else {
                bits_lsb_p1.insert(bits_lsb_p1.begin(), (assigns[i]==l_True ? '1' : '0'));
            }
        }

        // Checking if lsb bits of p2 equivalent to p_msb_count are set
        for(int i=q_lsb_var-1; i<q_lsb_var+p_msb_count-1; i++)
        {
            if(assigns[i]==l_Undef) {
                call_coppersmith_lsb_p2 = false;
                break;
            } else {
                bits_lsb_p2.insert(bits_lsb_p2.begin(), (assigns[i]==l_True ? '1' : '0'));
            }
        }

        if(call_coppersmith_lsb_p1 && opt_lo) {
            mpz_t mpz_p_tilda;
            mpz_init(mpz_p_tilda);

            for (char bit : bits_lsb_p1) {
                int bitValue = bit - '0';
                mpz_mul_ui(mpz_p_tilda, mpz_p_tilda, 2);
                mpz_add_ui(mpz_p_tilda, mpz_p_tilda, bitValue);
            }
            cs_count_lsb_p1++;
            double start = cpuTime();
            bool res = coppersmith(mpz_p_tilda, 1);
            double end = cpuTime();
            cs_time += end - start;
            if(res) {
                cout << "Coppersmith Successful!" << endl;
                cout << "Result from lo bits" << endl;
                cout << "Coppersmith was called " + to_string(cs_count) + " times." << endl;
                /*vec<Lit> clauses;
                size_t bitLength = mpz_sizeinbase(cs_prime, 2);
                int j = p_lsb_var-1;

                for (size_t i = 0; i < bitLength; ++i) {
                    clauses.clear();
                    int bitValue = mpz_tstbit(cs_prime, i);

                    clauses.push(mkLit(j++, bitValue?false:true));
                    addClause_(clauses);
                }*/
                printStats();
                printf("\nSATISFIABLE\n");
                exit(10);
            } 
            else {
                const int size = out_learnts.size();
                out_learnts.push();
                for(int i=p_lsb_var-1; i<p_lsb_var+p_msb_count-1; i++)
                {
                    out_learnts[size].push(mkLit(i, assigns[i]==l_True));
                }
            }
            mpz_clear(mpz_p_tilda);
        }

        if(call_coppersmith_lsb_p2 && opt_lo && opt_cs_both_primes) {
            mpz_t mpz_p_tilda;
            mpz_init(mpz_p_tilda);
            for (char bit : bits_lsb_p2) {
                int bitValue = bit - '0';
                mpz_mul_ui(mpz_p_tilda, mpz_p_tilda, 2);
                mpz_add_ui(mpz_p_tilda, mpz_p_tilda, bitValue);
            }
            cs_count_lsb_p2++;
            double start = cpuTime();
            bool res = coppersmith(mpz_p_tilda, 1);
            double end = cpuTime();
            cs_time += end - start;
            if(res) {
                cout << "Coppersmith Successful!" << endl;
                cout << "Result from lo bits" << endl;
                cout << "Coppersmith was called " + to_string(cs_count) + " times." << endl;
                /*vec<Lit> clauses;
                size_t bitLength = mpz_sizeinbase(cs_prime, 2);
                int j = q_lsb_var-1;

                for (size_t i = 0; i < bitLength; ++i) {
                    clauses.clear();
                    int bitValue = mpz_tstbit(cs_prime, i);

                    clauses.push(mkLit(j++, bitValue?false:true));
                    addClause_(clauses);
                }*/
                printStats();
                printf("\nSATISFIABLE\n");
                exit(10);
            } 
            else {
                const int size = out_learnts.size();
                out_learnts.push();
                for(int i=q_lsb_var-1; i<q_lsb_var+p_msb_count-1; i++)
                {
                    out_learnts[size].push(mkLit(i, assigns[i]==l_True));
                }
            }
            mpz_clear(mpz_p_tilda);
        }

        if(call_coppersmith_msb_p1 && opt_hi) {
            mpz_t mpz_p_tilda;
            mpz_init(mpz_p_tilda);
            for(int i=p_lsb_var-1; i<p_msb_var-p_msb_count; i++)
            {
                bits_msb_p1.push_back('0');
            }
            
            for (char bit : bits_msb_p1) {
                int bitValue = bit - '0';
                mpz_mul_ui(mpz_p_tilda, mpz_p_tilda, 2);
                mpz_add_ui(mpz_p_tilda, mpz_p_tilda, bitValue);
            }
            cs_count_msb_p1++;
            double start = cpuTime();
            bool res = coppersmith(mpz_p_tilda, 0);
            double end = cpuTime();
            cs_time += end - start;
            if(res) {
                cout << "Coppersmith Successful!" << endl;
                cout << "Result from hi bits" << endl;
                cout << "Coppersmith was called " + to_string(cs_count) + " times." << endl;
                printStats();
                printf("\nSATISFIABLE\n");
                exit(10);
            } 
            else {
                const int size = out_learnts.size();
                out_learnts.push();
                for(int i=p_msb_var-1; i>p_msb_var-p_msb_count-1; i--)
                {
                    out_learnts[size].push(mkLit(i, assigns[i]==l_True));
                }
            }
            mpz_clear(mpz_p_tilda);
        }

        if(call_coppersmith_msb_p2 && opt_hi && opt_cs_both_primes) {
            mpz_t mpz_p_tilda;
            mpz_init(mpz_p_tilda);
            for(int i=q_lsb_var-1; i<q_msb_var-p_msb_count; i++)
            {
                bits_msb_p2.push_back('0');
            }

            for (char bit : bits_msb_p2) {
                int bitValue = bit - '0';
                mpz_mul_ui(mpz_p_tilda, mpz_p_tilda, 2);
                mpz_add_ui(mpz_p_tilda, mpz_p_tilda, bitValue);
            }
            cs_count_msb_p2++;
            double start = cpuTime();
            bool res = coppersmith(mpz_p_tilda, 0);
            double end = cpuTime();
            cs_time += end - start;
            if(res) {
                cout << "Coppersmith Successful!" << endl;
                cout << "Result from hi bits" << endl;
                cout << "Coppersmith was called " + to_string(cs_count) + " times." << endl;
                printStats();
                printf("\nSATISFIABLE\n");
                exit(10);
            } 
            else {
                const int size = out_learnts.size();
                out_learnts.push();
                for(int i=q_msb_var-1; i>q_msb_var-p_msb_count-1; i--)
                {
                    out_learnts[size].push(mkLit(i, assigns[i]==l_True));
                }
            }
            mpz_clear(mpz_p_tilda);
        }
    }
    wait++;
}

bool Solver::assertingClause(CRef confl) {
    Clause& c = ca[confl];
    int asserting = -1;
    for (int i = 0; i < c.size(); i++) {
        if (value(c[i]) == l_Undef) {
            if (asserting != -1) return false;
            asserting = i;
        }
    }
    return asserting == 0;
}

void Solver::analyze(vec<Lit>& conflvec, vec<Lit>& out_learnt, int& out_btlevel)
{
    int pathC = 0;
    CRef confl;
    Lit p     = lit_Undef;

    int cur_max = level(var(conflvec[0]));
    for(int j=1; j < conflvec.size(); j++) {
        if(level(var(conflvec[j])) > cur_max) {
            cur_max = level(var(conflvec[j]));
        }
    }
    if(cur_max == 0) {
        out_btlevel = -1;
        return;
    }
    if (conflvec.size() == 1) {
        out_btlevel = 0;
        conflvec.copyTo(out_learnt);
        return;
    }

    // Generate conflict clause:
    //
    out_learnt.push();      // (leave room for the asserting literal)
    int index   = trail.size() - 1;

        for (int j = (p == lit_Undef) ? 0 : 1; j < conflvec.size(); j++){
            Lit q = conflvec[j];

            if (!seen[var(q)] && level(var(q)) > 0){
#if BRANCHING_HEURISTIC == CHB
                last_conflict[var(q)] = conflicts;
#elif BRANCHING_HEURISTIC == VSIDS
                varBumpActivity(var(q));
#endif
                conflicted[var(q)]++;
                seen[var(q)] = 1;
                if (level(var(q)) >= cur_max)
                    pathC++;
                else
                    out_learnt.push(q);
            }
        }

        // Select next clause to look at:
        while (!seen[var(trail[index--])]);
        p     = trail[index+1];
        confl = reason(var(p));
        seen[var(p)] = 0;
        pathC--;

    while(pathC > 0){
        assert(confl != CRef_Undef); // (otherwise should be UIP)
        Clause& c = ca[confl];

#if LBD_BASED_CLAUSE_DELETION
        if (c.learnt() && c.activity() > 2)
            c.activity() = lbd(c);
#else
        if (c.learnt())
            claBumpActivity(c);
#endif

        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
            Lit q = c[j];

            if (!seen[var(q)] && level(var(q)) > 0){
#if BRANCHING_HEURISTIC == CHB
                last_conflict[var(q)] = conflicts;
#elif BRANCHING_HEURISTIC == VSIDS
                varBumpActivity(var(q));
#endif
                conflicted[var(q)]++;
                seen[var(q)] = 1;
                if (level(var(q)) >= cur_max)
                    pathC++;
                else
                    out_learnt.push(q);
            }
        }

        // Select next clause to look at:
        while (!seen[var(trail[index--])]);
        p     = trail[index+1];
        confl = reason(var(p));
        seen[var(p)] = 0;
        pathC--;

    }
    out_learnt[0] = ~p;

    // Simplify conflict clause:
    //
    int i, j;
    out_learnt.copyTo(analyze_toclear);
    if (ccmin_mode == 2){
        uint32_t abstract_level = 0;
        for (i = 1; i < out_learnt.size(); i++)
            abstract_level |= abstractLevel(var(out_learnt[i])); // (maintain an abstraction of levels involved in conflict)

        for (i = j = 1; i < out_learnt.size(); i++)
            if (reason(var(out_learnt[i])) == CRef_Undef || !litRedundant(out_learnt[i], abstract_level))
                out_learnt[j++] = out_learnt[i];

    }else if (ccmin_mode == 1){
        for (i = j = 1; i < out_learnt.size(); i++){
            Var x = var(out_learnt[i]);

            if (reason(x) == CRef_Undef)
                out_learnt[j++] = out_learnt[i];
            else{
                Clause& c = ca[reason(var(out_learnt[i]))];
                for (int k = 1; k < c.size(); k++)
                    if (!seen[var(c[k])] && level(var(c[k])) > 0){
                        out_learnt[j++] = out_learnt[i];
                        break; }
            }
        }
    }else
        i = j = out_learnt.size();

    max_literals += out_learnt.size();
    out_learnt.shrink(i - j);
    tot_literals += out_learnt.size();

    // Find correct backtrack level:
    //
    if (out_learnt.size() == 1)
        out_btlevel = 0;
    else{
        int max_i = 1;
        // Find the first literal assigned at the next-highest level:
        for (int i = 2; i < out_learnt.size(); i++)
            if (level(var(out_learnt[i])) > level(var(out_learnt[max_i])))
                max_i = i;
        // Swap-in this literal at index 1:
        Lit p             = out_learnt[max_i];
        out_learnt[max_i] = out_learnt[1];
        out_learnt[1]     = p;
        out_btlevel       = level(var(p));
    }

#if ALMOST_CONFLICT
    seen[var(p)] = true;
    for(int i = out_learnt.size() - 1; i >= 0; i--) {
        Var v = var(out_learnt[i]);
        CRef rea = reason(v);
        if (rea != CRef_Undef) {
            Clause& reaC = ca[rea];
            for (int i = 0; i < reaC.size(); i++) {
                Lit l = reaC[i];
                if (!seen[var(l)]) {
                    seen[var(l)] = true;
                    almost_conflicted[var(l)]++;
                    analyze_toclear.push(l);
                }
            }
        }
    }
#endif
    for (int j = 0; j < analyze_toclear.size(); j++) seen[var(analyze_toclear[j])] = 0;    // ('seen[]' is now cleared)
}

/*_________________________________________________________________________________________________
|
|  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
|  
|  Description:
|    Analyze conflict and produce a reason clause.
|  
|    Pre-conditions:
|      * 'out_learnt' is assumed to be cleared.
|      * Current decision level must be greater than root level.
|  
|    Post-conditions:
|      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
|      * If out_learnt.size() > 1 then 'out_learnt[1]' has the greatest decision level of the 
|        rest of literals. There may be others from the same level though.
|  
|________________________________________________________________________________________________@*/
void Solver::analyze(CRef confl, vec<Lit>& out_learnt, int& out_btlevel)
{
    int pathC = 0;
    Lit p     = lit_Undef;

    // Generate conflict clause:
    //
    out_learnt.push();      // (leave room for the asserting literal)
    int index   = trail.size() - 1;

    do{
        assert(confl != CRef_Undef); // (otherwise should be UIP)
        Clause& c = ca[confl];

#if LBD_BASED_CLAUSE_DELETION
        if (c.learnt() && c.activity() > 2)
            c.activity() = lbd(c);
#else
        if (c.learnt())
            claBumpActivity(c);
#endif

        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
            Lit q = c[j];

            if (!seen[var(q)] && level(var(q)) > 0){
#if BRANCHING_HEURISTIC == CHB
                last_conflict[var(q)] = conflicts;
#elif BRANCHING_HEURISTIC == VSIDS
                varBumpActivity(var(q));
#endif
                conflicted[var(q)]++;
                seen[var(q)] = 1;
                if (level(var(q)) >= decisionLevel())
                    pathC++;
                else
                    out_learnt.push(q);
            }
        }
        
        // Select next clause to look at:
        while (!seen[var(trail[index--])]);
        p     = trail[index+1];
        confl = reason(var(p));
        seen[var(p)] = 0;
        pathC--;

    }while (pathC > 0);
    out_learnt[0] = ~p;

    // Simplify conflict clause:
    //
    int i, j;
    out_learnt.copyTo(analyze_toclear);
    if (ccmin_mode == 2){
        uint32_t abstract_level = 0;
        for (i = 1; i < out_learnt.size(); i++)
            abstract_level |= abstractLevel(var(out_learnt[i])); // (maintain an abstraction of levels involved in conflict)

        for (i = j = 1; i < out_learnt.size(); i++)
            if (reason(var(out_learnt[i])) == CRef_Undef || !litRedundant(out_learnt[i], abstract_level))
                out_learnt[j++] = out_learnt[i];
        
    }else if (ccmin_mode == 1){
        for (i = j = 1; i < out_learnt.size(); i++){
            Var x = var(out_learnt[i]);

            if (reason(x) == CRef_Undef)
                out_learnt[j++] = out_learnt[i];
            else{
                Clause& c = ca[reason(var(out_learnt[i]))];
                for (int k = 1; k < c.size(); k++)
                    if (!seen[var(c[k])] && level(var(c[k])) > 0){
                        out_learnt[j++] = out_learnt[i];
                        break; }
            }
        }
    }else
        i = j = out_learnt.size();

    max_literals += out_learnt.size();
    out_learnt.shrink(i - j);
    tot_literals += out_learnt.size();

    // Find correct backtrack level:
    //
    if (out_learnt.size() == 1)
        out_btlevel = 0;
    else{
        int max_i = 1;
        // Find the first literal assigned at the next-highest level:
        for (int i = 2; i < out_learnt.size(); i++)
            if (level(var(out_learnt[i])) > level(var(out_learnt[max_i])))
                max_i = i;
        // Swap-in this literal at index 1:
        Lit p             = out_learnt[max_i];
        out_learnt[max_i] = out_learnt[1];
        out_learnt[1]     = p;
        out_btlevel       = level(var(p));
    }

#if ALMOST_CONFLICT
    seen[var(p)] = true;
    for(int i = out_learnt.size() - 1; i >= 0; i--) {
        Var v = var(out_learnt[i]);
        CRef rea = reason(v);
        if (rea != CRef_Undef) {
            Clause& reaC = ca[rea];
            for (int i = 0; i < reaC.size(); i++) {
                Lit l = reaC[i];
                if (!seen[var(l)]) {
                    seen[var(l)] = true;
                    almost_conflicted[var(l)]++;
                    analyze_toclear.push(l);
                }
            }
        }
    }
#endif
    for (int j = 0; j < analyze_toclear.size(); j++) seen[var(analyze_toclear[j])] = 0;    // ('seen[]' is now cleared)
}


// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
bool Solver::litRedundant(Lit p, uint32_t abstract_levels)
{
    analyze_stack.clear(); analyze_stack.push(p);
    int top = analyze_toclear.size();
    while (analyze_stack.size() > 0){
        assert(reason(var(analyze_stack.last())) != CRef_Undef);
        Clause& c = ca[reason(var(analyze_stack.last()))]; analyze_stack.pop();

        for (int i = 1; i < c.size(); i++){
            Lit p  = c[i];
            if (!seen[var(p)] && level(var(p)) > 0){
                if (reason(var(p)) != CRef_Undef && (abstractLevel(var(p)) & abstract_levels) != 0){
                    seen[var(p)] = 1;
                    analyze_stack.push(p);
                    analyze_toclear.push(p);
                }else{
                    for (int j = top; j < analyze_toclear.size(); j++)
                        seen[var(analyze_toclear[j])] = 0;
                    analyze_toclear.shrink(analyze_toclear.size() - top);
                    return false;
                }
            }
        }
    }

    return true;
}


/*_________________________________________________________________________________________________
|
|  analyzeFinal : (p : Lit)  ->  [void]
|  
|  Description:
|    Specialized analysis procedure to express the final conflict in terms of assumptions.
|    Calculates the (possibly empty) set of assumptions that led to the assignment of 'p', and
|    stores the result in 'out_conflict'.
|________________________________________________________________________________________________@*/
void Solver::analyzeFinal(Lit p, vec<Lit>& out_conflict)
{
    out_conflict.clear();
    out_conflict.push(p);

    if (decisionLevel() == 0)
        return;

    seen[var(p)] = 1;

    for (int i = trail.size()-1; i >= trail_lim[0]; i--){
        Var x = var(trail[i]);
        if (seen[x]){
            if (reason(x) == CRef_Undef){
                assert(level(x) > 0);
                out_conflict.push(~trail[i]);
            }else{
                Clause& c = ca[reason(x)];
                for (int j = 1; j < c.size(); j++)
                    if (level(var(c[j])) > 0)
                        seen[var(c[j])] = 1;
            }
            seen[x] = 0;
        }
    }

    seen[var(p)] = 0;
}


void Solver::uncheckedEnqueue(Lit p, CRef from)
{
    assert(value(p) == l_Undef);
    picked[var(p)] = conflicts;
#if ANTI_EXPLORATION
    uint64_t age = conflicts - canceled[var(p)];
    if (age > 0) {
        double decay = pow(0.95, age);
        activity[var(p)] *= decay;
        if (order_heap.inHeap(var(p))) {
            order_heap.increase(var(p));
        }
    }
#endif
    conflicted[var(p)] = 0;
#if ALMOST_CONFLICT
    almost_conflicted[var(p)] = 0;
#endif
    assigns[var(p)] = lbool(!sign(p));
    vardata[var(p)] = mkVarData(from, decisionLevel());
    trail.push_(p);
}


/*_________________________________________________________________________________________________
|
|  propagate : [void]  ->  [Clause*]
|  
|  Description:
|    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
|    otherwise CRef_Undef.
|  
|    Post-conditions:
|      * the propagation queue is empty, even if there was a conflict.
|________________________________________________________________________________________________@*/
CRef Solver::propagate()
{
    CRef    confl     = CRef_Undef;
    int     num_props = 0;
    watches.cleanAll();

    while (qhead < trail.size()){
        Lit            p   = trail[qhead++];     // 'p' is enqueued fact to propagate.
        vec<Watcher>&  ws  = watches[p];
        Watcher        *i, *j, *end;
        num_props++;

        for (i = j = (Watcher*)ws, end = i + ws.size();  i != end;){
            // Try to avoid inspecting the clause:
            Lit blocker = i->blocker;
            if (value(blocker) == l_True){
                *j++ = *i++; continue; }

            // Make sure the false literal is data[1]:
            CRef     cr        = i->cref;
            Clause&  c         = ca[cr];
            Lit      false_lit = ~p;
            if (c[0] == false_lit)
                c[0] = c[1], c[1] = false_lit;
            assert(c[1] == false_lit);
            i++;

            // If 0th watch is true, then clause is already satisfied.
            Lit     first = c[0];
            Watcher w     = Watcher(cr, first);
            if (first != blocker && value(first) == l_True){
                *j++ = w; continue; }

            // Look for new watch:
            for (int k = 2; k < c.size(); k++)
                if (value(c[k]) != l_False){
                    c[1] = c[k]; c[k] = false_lit;
                    watches[~c[1]].push(w);
                    goto NextClause; }

            // Did not find watch -- clause is unit under assignment:
            *j++ = w;
            if (value(first) == l_False){
                confl = cr;
                qhead = trail.size();
                // Copy the remaining watches:
                while (i < end)
                    *j++ = *i++;
            }else
                uncheckedEnqueue(first, cr);

        NextClause:;
        }
        ws.shrink(i - j);
    }
    propagations += num_props;
    simpDB_props -= num_props;

    return confl;
}

int min(int a, int b) {
    return a < b ? a : b;
}

/*_________________________________________________________________________________________________
|
|  reduceDB : ()  ->  [void]
|  
|  Description:
|    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
|    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
|________________________________________________________________________________________________@*/
struct reduceDB_lt { 
    ClauseAllocator& ca;
#if LBD_BASED_CLAUSE_DELETION
    vec<double>& activity;
    reduceDB_lt(ClauseAllocator& ca_,vec<double>& activity_) : ca(ca_), activity(activity_) {}
#else
    reduceDB_lt(ClauseAllocator& ca_) : ca(ca_) {}
#endif
    bool operator () (CRef x, CRef y) { 
#if LBD_BASED_CLAUSE_DELETION
        return ca[x].activity() > ca[y].activity();
    }
#else
        return ca[x].size() > 2 && (ca[y].size() == 2 || ca[x].activity() < ca[y].activity()); } 
#endif
};
void Solver::reduceDB()
{
    int     i, j;
#if LBD_BASED_CLAUSE_DELETION
    sort(learnts, reduceDB_lt(ca, activity));
#else
    double  extra_lim = cla_inc / learnts.size();    // Remove any clause below this activity
    sort(learnts, reduceDB_lt(ca));
#endif

    // Don't delete binary or locked clauses. From the rest, delete clauses from the first half
    // and clauses with activity smaller than 'extra_lim':
#if LBD_BASED_CLAUSE_DELETION
    for (i = j = 0; i < learnts.size(); i++){
        Clause& c = ca[learnts[i]];
        if (c.activity() > 2 && !locked(c) && i < learnts.size() / 2)
#else
    for (i = j = 0; i < learnts.size(); i++){
        Clause& c = ca[learnts[i]];
        if (c.size() > 2 && !locked(c) && (i < learnts.size() / 2 || c.activity() < extra_lim))
#endif
            removeClause(learnts[i]);
        else
            learnts[j++] = learnts[i];
    }
    learnts.shrink(i - j);
    checkGarbage();
}


void Solver::removeSatisfied(vec<CRef>& cs)
{
    int i, j;
    for (i = j = 0; i < cs.size(); i++){
        Clause& c = ca[cs[i]];
        if (satisfied(c))
            removeClause(cs[i]);
        else
            cs[j++] = cs[i];
    }
    cs.shrink(i - j);
}


void Solver::rebuildOrderHeap()
{
    vec<Var> vs;
    for (Var v = 0; v < nVars(); v++)
        if (decision[v] && value(v) == l_Undef)
            vs.push(v);
    order_heap.build(vs);
}


/*_________________________________________________________________________________________________
|
|  simplify : [void]  ->  [bool]
|  
|  Description:
|    Simplify the clause database according to the current top-level assigment. Currently, the only
|    thing done here is the removal of satisfied clauses, but more things can be put here.
|________________________________________________________________________________________________@*/
bool Solver::simplify()
{
    assert(decisionLevel() == 0);

    if (!ok || propagate() != CRef_Undef)
        return ok = false;

    if (nAssigns() == simpDB_assigns || (simpDB_props > 0))
        return true;

    // Remove satisfied clauses:
    removeSatisfied(learnts);
    if (remove_satisfied)        // Can be turned off.
        removeSatisfied(clauses);
    checkGarbage();
    rebuildOrderHeap();

    simpDB_assigns = nAssigns();
    simpDB_props   = clauses_literals + learnts_literals;   // (shouldn't depend on stats really, but it will do for now)

    return true;
}

/*_________________________________________________________________________________________________
|
|  search : (nof_conflicts : int) (params : const SearchParams&)  ->  [lbool]
|  
|  Description:
|    Search for a model the specified number of conflicts. 
|    NOTE! Use negative value for 'nof_conflicts' indicate infinity.
|  
|  Output:
|    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
|    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
|    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
|________________________________________________________________________________________________@*/
lbool Solver::search(int nof_conflicts)
{
    assert(ok);
    int         backtrack_level;
    int         conflictC = 0;
    vec<Lit>    learnt_clause;
    vec<Lit>    units;
    starts++;

    for (;;){
        CRef confl = propagate();

#if BRANCHING_HEURISTIC == CHB
        double multiplier = confl == CRef_Undef ? reward_multiplier : 1.0;
        for (int a = action; a < trail.size(); a++) {
            Var v = var(trail[a]);
            uint64_t age = conflicts - last_conflict[v] + 1;
            double reward = multiplier / age ;
            double old_activity = activity[v];
            activity[v] = step_size * reward + ((1 - step_size) * old_activity);
            if (order_heap.inHeap(v)) {
                if (activity[v] > old_activity)
                    order_heap.decrease(v);
                else
                    order_heap.increase(v);
            }
        }
#endif
        if (confl != CRef_Undef){
            // CONFLICT
            conflicts++; conflictC++;
#if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB
            if (step_size > min_step_size)
                step_size -= step_size_dec;
#endif
            if (decisionLevel() == 0) return l_False;

            learnt_clause.clear();
            analyze(confl, learnt_clause, backtrack_level);

            cancelUntil(backtrack_level);

#if BRANCHING_HEURISTIC == CHB
            action = trail.size();
#endif

            if (learnt_clause.size() == 1){
                uncheckedEnqueue(learnt_clause[0]);
            }else{
                CRef cr = ca.alloc(learnt_clause, true);
                learnts.push(cr);
                attachClause(cr);
#if LBD_BASED_CLAUSE_DELETION
                Clause& clause = ca[cr];
                clause.activity() = lbd(clause);
#else
                claBumpActivity(ca[cr]);
#endif
                uncheckedEnqueue(learnt_clause[0], cr);
            }

#if BRANCHING_HEURISTIC == VSIDS
            varDecayActivity();
#endif
#if ! LBD_BASED_CLAUSE_DELETION
            claDecayActivity();
#endif

            if (--learntsize_adjust_cnt == 0){
                learntsize_adjust_confl *= learntsize_adjust_inc;
                learntsize_adjust_cnt    = (int)learntsize_adjust_confl;
#if ! RAPID_DELETION
                max_learnts             *= learntsize_inc;
#endif

                if (verbosity >= 1) {
                    printf("| %9d | %7d %8d %8d | %8d %8d %6.0f | %6.3f %% |",
                               (int)conflicts, 
                               (int)dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]), nClauses(), (int)clauses_literals, 
                               (int)max_learnts, nLearnts(), (double)learnts_literals/nLearnts(), progressEstimate()*100);
                    if(!opt_only_sat && (opt_lo || opt_hi)) printf(" %10ld |", cs_count);
                    if(!opt_only_sat && opt_prog_hc) printf(" %10ld |", hc_count);
                    printf("\n");
                    fflush(stdout);
                }
            }

        }else{
            // NO CONFLICT
            if (nof_conflicts >= 0 && conflictC >= nof_conflicts || !withinBudget()){
                // Reached bound on number of conflicts:
                progress_estimate = progressEstimate();
                cancelUntil(0);
                return l_Undef; }

            // Simplify the set of problem clauses:
            if (decisionLevel() == 0 && !simplify())
                return l_False;

            if (learnts.size()-nAssigns() >= max_learnts) {
                // Reduce the set of learnt clauses:
                reduceDB();
#if RAPID_DELETION
                max_learnts += 500;
#endif
            }

            Lit next = lit_Undef;
            Lit branch_lit = lit_Undef;
            while (decisionLevel() < assumptions.size()){
                // Perform user provided assumption:
                Lit p = assumptions[decisionLevel()];
                if (value(p) == l_True){
                    // Dummy decision level:
                    newDecisionLevel();
                }else if (value(p) == l_False){
                    analyzeFinal(~p, conflict);
                    return l_False;
                }else{
                    next = p;
                    break;
                }
            }

            if (next == lit_Undef){
                // New variable decision:
                decisions++;
                next = pickBranchLit();

                if(verbosity >= 3) {
                    // Print the bits of the primes if verb=2
                    for(int i=prime_len; i>0; i--) {
                        printf("%c", assigns[i-1]==l_True ? '1' : assigns[i-1]==l_False ? '0' : '?');
                    }
                    printf(" ");
                    for(int i=q_msb_var; i>q_msb_var-prime_len; i--) {
                        printf("%c", assigns[i-1]==l_True ? '1' : assigns[i-1]==l_False ? '0' : '?');
                    }
                    printf("\n");
                }

                callbackLearntClauses.clear();
                callbackFunction(next == lit_Undef, callbackLearntClauses, branch_lit);
                if (callbackLearntClauses.size() > 0) {
                    conflicts++; conflictC++;
                    int pending = learnts.size();
                    units.clear();
                    backtrack_level = decisionLevel();
                    if(verbosity >= 2)
                        printf("%d learned clauses\n", callbackLearntClauses.size());
                    for (int i = 0; i < callbackLearntClauses.size(); i++) {
                        if(verbosity >= 2)
                            printClause(callbackLearntClauses[i]);
                        int level;
                        learnt_clause.clear();
                        analyze(callbackLearntClauses[i], learnt_clause, level);
                        if (level == -1) {
                            return l_False;
                        } else if (level < backtrack_level) {
                            backtrack_level = level;
                        }
                        if (learnt_clause.size() == 1) {
                            units.push(learnt_clause[0]);
                        } else {
                            CRef cr = ca.alloc(learnt_clause, true);
                            learnts.push(cr);
                            attachClause(cr);
#if LBD_BASED_CLAUSE_DELETION
                            Clause& clause = ca[cr];
                            clause.activity() = lbd(clause);
#else
                            claBumpActivity(ca[cr]);
#endif
                        }
                    }

                    cancelUntil(backtrack_level);

#if BRANCHING_HEURISTIC == CHB
                    action = trail.size();
#endif

                    for (int i = 0; i < units.size(); i++) {
                        Lit l = units[i];
                        // Make sure it wasn't assigned by one of the other callback learnt clauses.
                        if (value(l) == l_Undef) uncheckedEnqueue(l);
                    }
                    for (int i = pending; i < learnts.size(); i++) {
                        CRef cr = learnts[i];
                        Clause& c = ca[cr];
                        bool asserting = assertingClause(cr);
                        if (asserting) uncheckedEnqueue(c[0], cr);
                    }
                    // Do not branch.
                    if (next != lit_Undef) {
                        insertVarOrder(var(next));
                        next = lit_Undef;
                    }
                } else if (next == lit_Undef)
                    // Model found:
                    return l_True;
            }

            if (branch_lit != lit_Undef) {
                // if(value(branch_lit) != l_Undef) printf("Error: branching variable must be unassigned\n");
                // printf("Branching on literal %s%d\n", sign(branch_lit) ? "-" : "", var(branch_lit)+1);
                // Increase decision level and enqueue 'branch_lit'
                newDecisionLevel();
    #if BRANCHING_HEURISTIC == CHB
                action = trail.size();
    #endif
                uncheckedEnqueue(branch_lit);
            } else if (next != lit_Undef) {
                // Increase decision level and enqueue 'next'
                newDecisionLevel();
    #if BRANCHING_HEURISTIC == CHB
                action = trail.size();
    #endif
                uncheckedEnqueue(next);
            }
        }
    }
}


double Solver::progressEstimate() const
{
    double  progress = 0;
    double  F = 1.0 / nVars();

    for (int i = 0; i <= decisionLevel(); i++){
        int beg = i == 0 ? 0 : trail_lim[i - 1];
        int end = i == decisionLevel() ? trail.size() : trail_lim[i];
        progress += pow(F, i) * (end - beg);
    }

    return progress / nVars();
}

/*
  Finite subsequences of the Luby-sequence:

  0: 1
  1: 1 1 2
  2: 1 1 2 1 1 2 4
  3: 1 1 2 1 1 2 4 1 1 2 1 1 2 4 8
  ...


 */

static double luby(double y, int x){

    // Find the finite subsequence that contains index 'x', and the
    // size of that subsequence:
    int size, seq;
    for (size = 1, seq = 0; size < x+1; seq++, size = 2*size+1);

    while (size-1 != x){
        size = (size-1)>>1;
        seq--;
        x = x % size;
    }

    return pow(y, seq);
}

// NOTE: assumptions passed in member-variable 'assumptions'.
lbool Solver::solve_()
{
    model.clear();
    conflict.clear();
    if (!ok) return l_False;

    solves++;

#if RAPID_DELETION
    max_learnts               = 2000;
#else
    max_learnts               = nClauses() * learntsize_factor;
#endif
    learntsize_adjust_confl   = learntsize_adjust_start_confl;
    learntsize_adjust_cnt     = (int)learntsize_adjust_confl;
    lbool   status            = l_Undef;

    if (verbosity >= 1){
        if (verbosity >= 2){
            printf("LBD Based Clause Deletion : %d\n", LBD_BASED_CLAUSE_DELETION);
            printf("Rapid Deletion : %d\n", RAPID_DELETION);
            printf("Almost Conflict : %d\n", ALMOST_CONFLICT);
            printf("Anti Exploration : %d\n", ANTI_EXPLORATION);
        }
        printf("============================[ Search Statistics ]=============================="); if(!opt_only_sat) { if(opt_lo || opt_hi) printf("============="); if(opt_prog_hc) printf("============="); } printf("\n");
        printf("| Conflicts |          ORIGINAL         |          LEARNT          | Progress |"); if(!opt_only_sat) { if(opt_lo || opt_hi) printf(" CS Clauses |"); if(opt_prog_hc) printf(" HC Branch  |"); } printf("\n");
        printf("|           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |          |"); if(!opt_only_sat) { if(opt_lo || opt_hi) printf("            |"); if(opt_prog_hc) printf("            |"); } printf("\n");
        printf("==============================================================================="); if(!opt_only_sat) { if(opt_lo || opt_hi) printf("============="); if(opt_prog_hc) printf("============="); } printf("\n");
    }

    // Search:
    int curr_restarts = 0;
    while (status == l_Undef){
        double rest_base = luby_restart ? luby(restart_inc, curr_restarts) : pow(restart_inc, curr_restarts);
        status = search(rest_base * restart_first);
        if (!withinBudget()) break;
        curr_restarts++;
    }

    if (verbosity >= 1) {
        printf("==============================================================================="); if(!opt_only_sat) { if(opt_lo || opt_hi) printf("============="); if(opt_prog_hc) printf("============="); } printf("\n");
    }

    if (status == l_True){
        // Extend & copy model:
        model.growTo(nVars());
        for (int i = 0; i < nVars(); i++) model[i] = value(i);
    }else if (status == l_False && conflict.size() == 0)
        ok = false;

    cancelUntil(0);
    return status;
}

//=================================================================================================
// Writing CNF to DIMACS:
// 
// FIXME: this needs to be rewritten completely.

static Var mapVar(Var x, vec<Var>& map, Var& max)
{
    if (map.size() <= x || map[x] == -1){
        map.growTo(x+1, -1);
        map[x] = max++;
    }
    return map[x];
}


void Solver::toDimacs(FILE* f, Clause& c, vec<Var>& map, Var& max)
{
    if (satisfied(c)) return;

    for (int i = 0; i < c.size(); i++)
        if (value(c[i]) != l_False)
            fprintf(f, "%s%d ", sign(c[i]) ? "-" : "", mapVar(var(c[i]), map, max)+1);
    fprintf(f, "0\n");
}


void Solver::toDimacs(const char *file, const vec<Lit>& assumps)
{
    FILE* f = fopen(file, "wr");
    if (f == NULL)
        fprintf(stderr, "could not open file %s\n", file), exit(1);
    toDimacs(f, assumps);
    fclose(f);
}


void Solver::toDimacs(FILE* f, const vec<Lit>& assumps)
{
    // Handle case when solver is in contradictory state:
    if (!ok){
        fprintf(f, "p cnf 1 2\n1 0\n-1 0\n");
        return; }

    vec<Var> map; Var max = 0;

    // Cannot use removeClauses here because it is not safe
    // to deallocate them at this point. Could be improved.
    int cnt = 0;
    for (int i = 0; i < clauses.size(); i++)
        if (!satisfied(ca[clauses[i]]))
            cnt++;
        
    for (int i = 0; i < clauses.size(); i++)
        if (!satisfied(ca[clauses[i]])){
            Clause& c = ca[clauses[i]];
            for (int j = 0; j < c.size(); j++)
                if (value(c[j]) != l_False)
                    mapVar(var(c[j]), map, max);
        }

    // Assumptions are added as unit clauses:
    cnt += assumptions.size();

    fprintf(f, "p cnf %d %d\n", max, cnt);

    for (int i = 0; i < assumptions.size(); i++){
        assert(value(assumptions[i]) != l_False);
        fprintf(f, "%s%d 0\n", sign(assumptions[i]) ? "-" : "", mapVar(var(assumptions[i]), map, max)+1);
    }

    for (int i = 0; i < clauses.size(); i++)
        toDimacs(f, ca[clauses[i]], map, max);

    if (verbosity > 0)
        printf("Wrote %d clauses with %d variables.\n", cnt, max);
}


//=================================================================================================
// Garbage Collection methods:

void Solver::relocAll(ClauseAllocator& to)
{
    // All watchers:
    //
    // for (int i = 0; i < watches.size(); i++)
    watches.cleanAll();
    for (int v = 0; v < nVars(); v++)
        for (int s = 0; s < 2; s++){
            Lit p = mkLit(v, s);
            // printf(" >>> RELOCING: %s%d\n", sign(p)?"-":"", var(p)+1);
            vec<Watcher>& ws = watches[p];
            for (int j = 0; j < ws.size(); j++)
                ca.reloc(ws[j].cref, to);
        }

    // All reasons:
    //
    for (int i = 0; i < trail.size(); i++){
        Var v = var(trail[i]);

        if (reason(v) != CRef_Undef && (ca[reason(v)].reloced() || locked(ca[reason(v)])))
            ca.reloc(vardata[v].reason, to);
    }

    // All learnt:
    //
    for (int i = 0; i < learnts.size(); i++)
        ca.reloc(learnts[i], to);

    // All original:
    //
    for (int i = 0; i < clauses.size(); i++)
        ca.reloc(clauses[i], to);
}

void Solver::garbageCollect()
{
    // Initialize the next region to a size corresponding to the estimated utilization degree. This
    // is not precise but should avoid some unnecessary reallocations for the new region:
    ClauseAllocator to(ca.size() - ca.wasted());
    relocAll(to);
    if (verbosity >= 2)
        printf("|  Garbage collection:   %12d bytes => %12d bytes             |\n", 
               ca.size()*ClauseAllocator::Unit_Size, to.size()*ClauseAllocator::Unit_Size);
    to.moveTo(ca);
}
