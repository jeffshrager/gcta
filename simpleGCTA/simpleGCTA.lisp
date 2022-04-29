; (load (compile-file "simpleGCTA"))

#|

Health is a real value between 1 (healthy) and 0 (dead).

As with the GCTA model, here each patient and each drug has a binary
feature vector associated with it, and the effectiveness of a
treatment depends upon the match between these. If you treat a
non-mutation you get a negative effect. For the moment, the effect
sizes, which are per-cycle (e.g., per month?) are the same.

Cohorts are groups of patients being treated the same. HOWEVER,
cohorts also are defined by a pattern describing patients in this
cohort in terms of a trinary f.v., where we also have ? for a don't
care (1 or 0) in that position. So you can have patients with
different mutation patterns in the same cohort, by virtue of the
?s. Moreover, you can have different cohorts with the same pattern,
but being treated with a different drug. Cohorts can also be defined
by treatement history, which is given by another 3ary f.v., which is
this time applied to the drugs that patients can have been exposed to
in this cohort -- 1=must have had this drug, 0=must not have had this
drug, ?=don't care. (It may be that we can use a pt-fv-length'ed 3ary
f.v. for the drugs. (Note that there is a subtle ... difference?
semi-intentional confusion? ... between a drug and a drug
target. There actually aren't drugs, per se, there are only targets,
and we are assuming that the system knows (in the omnipotent sense
of "know") which targets a given tx has (by virtue of the tx's
f.v.). So, when a pt has a history of particular txs, what they really
have is the union of the targets that various txs that they have had
contiains.xx

With all these variables, the number of cohorts can be extremely
large. Searching this space efficiently is the whole goal!

For the moment we don't deal with combinations of drugs. (Actually,
beacuse the drug f.v. can have multiple 1s, it is almost like a drug
whose f.v. has more than one 1 is actually a cocktail, but there are
other things that a cocktail may offer or disallow, which this
simplification does not capture.)

|#

