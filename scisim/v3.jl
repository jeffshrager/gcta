using Random, Distributions, Printf, StatsBase, Distances, Dates, SciPy
using PyPlot; const plt = PyPlot

# Exploring the iterated model. (See discussion in email from around
# Nov. 14th.) We start with a non-informative prior, and at each cycle
# there are a number of "experts" to choose from, including one
# special expert, who is the recipient (the external oncologist who is
# working with the patient). We have a model for each of these experts
# (which should eventually be learned from priors and data!) and we
# keep iterating through the expert assumed to be the best(for this
# case. Each time through we ask whether passing the report to the
# recipient (and it being processed through their expertise) will
# enable the case to reach the target "correctness threhsold" --
# nb. that is AFTER the recipient has thought about the case, so it's
# based on the same model as any expert, but the recipient is special
# because, even though s/he may not be the "best" expert on this case,
# we always prefer them if it will get us immediately below the
# correctness threhsold. Note that often we'll go way below the
# correctness threhsold. This seemingly anomalous behavior is because
# the recipient might not have got us to the threshold, so we passed
# it through an additional expert, which might have got us there, or
# very close, and THEN it's always passed through the recipient, often
# ended up WAY below the treshold.

# Both experts and cases are modeled by a 10 element binary feature
# vector. An expert's "expertise" for a given case is just the jaccard
# distance between her fvec and that of the case. Therefore, a better
# expert for a given case has a samller expertise value. When an
# expert is selected, that expertise is simple simply multiplied into
# the case score (which always starts at 1.0). 

global fvec_length = 10 # Applies to both experts and cases

mutable struct expert fvec end
new_expert()=expert([bitrand()[1] for i in 1:fvec_length])

# We create a bunch of experts globally, and then compute their
# distance to a given case at case consideration time.

global n_all_experts=100
global all_experts = [new_expert() for i in 1:n_all_experts]

struct case fvec end # Binary case feature vector (comes from ZachO via Bryan)
new_case()=case([bitrand()[1] for i in 1:fvec_length])

global target_score = 0.1

function run()
    global target_score, all_experts
    case = new_case()
    case_score = 1.0 # FFF tx_bv(...)
    expertises = [jaccard(case.fvec,expert.fvec) for expert in all_experts]
    recipient_expertise = jaccard(case.fvec,new_expert().fvec)
    while expected_tx_rank_score_after_expert_analysis(case_score,recipient_expertise)>=target_score
        (best_expected_score_post_consideration,best_expert_n) =
            findmin([expected_tx_rank_score_after_expert_analysis(case_score,expertise) for expertise in expertises])
        # UUU Once an expert is chosen, stop them from being chosen again by setting them to 1.0 UUU FFF !!!
        case_score = best_expected_score_post_consideration
        @printf("Expert #%s (@%s) brings the score to %s\n", best_expert_n,expertises[best_expert_n],case_score)
        expertises[best_expert_n]=1.0 # ??? What happens if we run out of experts?
    end
    final_score = round(recipient_expertise * case_score;sigdigits=2)
    # FFF Update the expert models based on various TBD analyses (esp. final recipient choice)
    @printf("Score released to recipent (with expertise: %s) bringing the final score to %s\n",
            round(recipient_expertise;sigdigits=2),final_score)
end

function expected_tx_rank_score_after_expert_analysis(case_score,expertise)
    return round(case_score * expertise;sigdigits=2)
end

# What we call a "treatment model" is a vector of normalized treatment
# features. This can be turned deterministically into a treatment
# ranking, however, experts traffic in the treatment model because it
# has features that are descriptive of multiple treatments (e.g.,
# "PD-1 inhibitor"), and the hypothesis is that experts reason at this
# level about treatments, not so much in the specific treatments
# (although specific treatments are part of the treatment model as
# well.

# We want to have a model of an expert which, given a particular
# ranking of treatments and a case (i.e., a feature vector for the
# case) tells us the expected improvement in ranking 

run()
