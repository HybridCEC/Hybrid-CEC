#ifndef _bump_h_INCLUDED
#define _bump_h_INCLUDED

struct kissat;

void kissat_bump_variables (struct kissat *);
void pol_dec_act(struct kissat *);
void bump_pol_sc (struct kissat *, unsigned);
void kissat_bump_one(struct kissat *);

#define MAX_SCORE 1e150

#endif
