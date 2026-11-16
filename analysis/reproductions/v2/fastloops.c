/* C inner loops for the repro v2 group A runs (P23, P27, P18).
 *
 * Each function performs the same arithmetic, in the same order, as the Python
 * model code it replaces; all random numbers are drawn in Python with the
 * model's own Generator. Bit-identity against the models is checked by
 * check_bit_identical.py. Compile without -ffast-math.
 */
#include <stdint.h>

/* --- Minority game (epc.models.minority_game.MinorityGame.simulate, S = 2) ---
 * scores   : (N, 2) float64, updated in place
 * strat_t  : (P, N, 2) int8, strat_t[mu, i, s] = strategies[i, s, mu]
 * coins    : (n_rounds, N) int64 coin tosses (used only if tie_random)
 * *mu_io   : history index on entry, updated on exit
 * att, win, mus : per-round outputs (attendance n1, winning side, mu used)
 */
void mg_rounds(double *scores, const int8_t *strat_t, int64_t N, int64_t P,
               const int64_t *coins, int tie_random, int linear, int64_t n_rounds,
               int64_t *mu_io, int64_t *att, int8_t *win, int64_t *mus)
{
    int64_t mu = *mu_io;
    int64_t half = N / 2;
    for (int64_t t = 0; t < n_rounds; t++) {
        const int8_t *sm = strat_t + mu * N * 2;
        int64_t n1 = 0;
        for (int64_t i = 0; i < N; i++) {
            double s0 = scores[2 * i], s1 = scores[2 * i + 1];
            int64_t best = (s1 > s0) ? 1 : 0;          /* np.argmax over two columns */
            if (tie_random && s0 == s1) best = coins[t * N + i];
            n1 += sm[2 * i + best];
        }
        int w = (n1 > half) ? 0 : 1;
        int64_t A = 2 * n1 - N;
        double f = linear ? (double)A : (A > 0 ? 1.0 : -1.0);
        for (int64_t i = 0; i < 2 * N; i++) {
            double a = 2.0 * (double)sm[i] - 1.0;      /* +-1 */
            scores[i] -= a * f;                        /* U_s += -a_s * (A or sign A) */
        }
        att[t] = n1;
        win[t] = (int8_t)w;
        mus[t] = mu;
        mu = ((mu << 1) | w) & (P - 1);
    }
    *mu_io = mu;
}

/* --- Nowak-May (epc.models.nowak_may.NowakMayModel.step, periodic) ---
 * grid, out : (R, C) int8 (0 = C, 1 = D); pay : (R, C) float64 scratch.
 * Returns the number of cooperators in out.
 */
int64_t nm_generation(const int8_t *grid, int8_t *out, double *pay,
                      int64_t R, int64_t C, double b)
{
    for (int64_t i = 0; i < R; i++) {
        for (int64_t j = 0; j < C; j++) {
            double cn = 0.0;
            for (int64_t dr = -1; dr <= 1; dr++) {
                int64_t r = i + dr; if (r < 0) r += R; else if (r >= R) r -= R;
                for (int64_t dc = -1; dc <= 1; dc++) {
                    int64_t c = j + dc; if (c < 0) c += C; else if (c >= C) c -= C;
                    cn += (double)(1 - grid[r * C + c]);
                }
            }
            pay[i * C + j] = (grid[i * C + j] == 0) ? cn * 1.0 : cn * b;
        }
    }
    int64_t ncoop = 0;
    for (int64_t i = 0; i < R; i++) {
        for (int64_t j = 0; j < C; j++) {
            double best = pay[i * C + j];
            int8_t bs = grid[i * C + j];
            for (int64_t dr = -1; dr <= 1; dr++) {
                int64_t r = i + dr; if (r < 0) r += R; else if (r >= R) r -= R;
                for (int64_t dc = -1; dc <= 1; dc++) {
                    if (dr == 0 && dc == 0) continue;
                    int64_t c = j + dc; if (c < 0) c += C; else if (c >= C) c -= C;
                    double p = pay[r * C + c];
                    if (p > best) { best = p; bs = grid[r * C + c]; }   /* strict > */
                }
            }
            out[i * C + j] = bs;
            ncoop += (bs == 0);
        }
    }
    return ncoop;
}

/* --- Voter (epc.models.voter.VoterModel.step) ---
 * Replaces: for k in range(N): grid[r, c] = grid[(r+dr) % R, (c+dc) % C]
 */
void voter_sweep(int8_t *grid, int64_t R, int64_t C,
                 const int64_t *rs, const int64_t *cs, const int64_t *nb,
                 int64_t n, const int64_t *dr, const int64_t *dc)
{
    for (int64_t k = 0; k < n; k++) {
        int64_t r = rs[k], c = cs[k], j = nb[k];
        int64_t nr = r + dr[j];
        int64_t nc = c + dc[j];
        if (nr < 0) nr += R; else if (nr >= R) nr -= R;
        if (nc < 0) nc += C; else if (nc >= C) nc -= C;
        grid[r * C + c] = grid[nr * C + nc];
    }
}
