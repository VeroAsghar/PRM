
using PRM

simp = SimParam(Tm=10, epsilon=1e-6, x0=0, h0=0.2, t0=0, tf=20, ts=1);

run_sim(simp)
