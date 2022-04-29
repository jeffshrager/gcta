using Random, Distributions, Printf, StatsBase, Distances, Dates
using PyPlot; const plt = PyPlot

# FFF Move items to a "bad" category
# %%% Remove the annoying search-by-sorting in choose_vtb_experts

#   Cases are a vector of 5 values in 0...9. This is also the "correct" treatment ranking for this case.
#   Experts specialities are also represented as a 5-vector in 0...9.
# (We use 0-9 instead of 1-10 so we can do compressed output of these vectors)
#   A recommendation (from a single expert) is a ranked 5-vector of treatments, also in 1-10, but
#     each also has a self-reported certainty. The closer the expert's specialization vector to
#     the patient's feature vector, the less "twiddled" the expert's opinion is (starting from the
#     true treatment vector, which is the same as the case vector -- this is probably too incestuous!)
#   In addition to the the expert's choices being twiddled more as they are farther from the patient's
#     feature vector, the less they self-weight their rankings. (That is, they provide an indication of
#     their overall certainty in their ranking which is essentially the same as the distance between the
#     patient's features/treatment vector and the expert's specializatiton vector.)
#   For each case, a VTB of 2 to 5 experts is convened, and they provide their rankings.

mutable struct case; features; end # Case is used elsewhere
mutable struct expert; name; spx; end # spx = specialization
mutable struct vtb; name; case; opinions; the_algorithm; algorithm_target_certainty; algorithm_recs;
                    inter_opinion_distances; algorithm_expert_distances; overall_score; end
mutable struct vtx; name; features; end

# ==================================================================
# These params need to be empirically adjusted to get the desired
# spread of expertises and treatments. ??? It's possible that the
# treatments and experts should be separately controlled.
# ==================================================================

# These control feature sampling for experts, treatments, and
# cases. You need to adjust them so that you get features centered on
# 5 with a reasonable range (e.g., the distribution report should look
# something like this: 0:6, 1:9, 2:40, 3:61, 4:90, 5:102, 6:83, 7:56,
# 8:34, 9:19) but also experts with a reasonable range of exterise
# vs. cases. (This is tested against treatments at the same time as
# the distrubution report is produced.) For this, a mean of 0.5 and a
# std of 0.2 is pretty good. (FFF It'd be nice to get the mean a bit
# higher, but then you wouldn't have a normal that tops at 1.0. As is,
# we don't have to use the cutoffs too much.) [??? !!! This isn't
# quite the metric we're interested in. We want something like a good
# distribution for each problem.]

global g_max_features_distance = 20.0
global g_feature_kurtosis = 2.0

# A high value of g_nexperts will slow things down, but will also
# increase the range of experts we can select from. Unfortunately that
# also makes it hard NOT to find a very good expert in any given case
# (all the cases will end up with a perfect expert). We do want
# experts to be good, but not always perfect! So you don't want to
# make this too high, or else you won't have a good range of
# expertise, but you also don't want it too low, or else you won't
# find good experts for very many cases.
global g_nexperts = 25 

# g_ntxs, the number of available treatments, doesn't seem to have
# much effect, based upon informal sensitivity analysis.
global g_ntxs = 100

# ==================================================================
# Utility stuff
# Randomizer initializations
global g_rng = MersenneTwister(1234);
global g_txindeces = collect(1:g_ntxs) # For Poisson sampling

global g_tracelevel = 0 
trace!(msg,lvl=0) = lvl >= g_tracelevel ? print(msg) : false

# ==================================================================
# The algorithm is one of the experts who has a skill level that is
# matched to the case at hand based on the level of learning that we
# believe we've undertaken (as opposed to based on
# specialization). ??? This has the (possible mis-)feature that the
# best the algorithm can ever be is as good as the best
# expert. There's other ways to do this. The best alternative is not
# to have the VTB have specializations at all, but rather to just use
# the algorithm_target_certainty directly to select algorithmic
# recommendations.
# ==================================================================

global g_experts = []
global g_txs = []
 
# ==================================================================
# Cases, expert's specialization, and treatments all use a 5-digit
# number. As each digit ranges 0-9 (in the usual way), and we center
# on 55555 (using randn and adding 5), so we're slightly skewed to the
# top. The parameters should be set so that we get mostly 5s, less 6s,
# and so on, and a few 9s. At the moment I sort of eyeball this...FFF
# this probably ought to be done in a mode principled way.
# ==================================================================

sample_features() = [max(0,(min(9,(convert(Int64,round(5+(randn(g_rng)*g_feature_kurtosis))))))) for i in 1:5]

a_case() = case(sample_features())
an_expert(name) = expert(name,sample_features())
a_tx(name) = vtx(name,sample_features())

function init_vtbs()
    global g_experts, g_nexperts, g_ntxs, g_txs
    g_experts = [an_expert(i) for i in 1:g_nexperts]
    g_txs = [a_tx(i) for i in 1:g_ntxs]
    @printf("Treatement feature distributions:\n")
    for i in 0:9; @printf("%s:%s, ",i,sum([count(e->(e==i),tx.features) for tx in g_txs])); end
    print("\n")
    # ??? !!! This isn't quite the metric we're interested in. We want
    # something like a good distribution for each problem.
    exvtx = [expert_level(expert,tx) for expert in g_experts for tx in g_txs]
    @printf("Experts vs. Cases (using treatments as a stand in for cases): mean=%s, std=%s\n", round(mean(exvtx);sigdigits=2),round(std(exvtx);sigdigits=2))
end

# ==================================================================
# The way we calculate distance to a treatment is by counting up the
# distance that each is apart from the other along each dimension, so,
# for example, [1,2,3,4,5] vs. [9,8,7,6,5] has distances:
# [9-1=8,...5-5=0] and we add all those up (their abs's, of course!)
# Unfortunately, with this method the range of possible exertise is so
# large that we need thousands of experts in order to find ones that
# are any good at all. Therefore, what we do is to use a simple
# gaussian sampling, centered at 5, to choose the
# ==================================================================

feature_distance(a,b) = sum([abs(a[i]-b[i]) for i in 1:5])

# ==================================================================
# Note that "expert_level" is always relative to the case at hand;
# It's NOT a property of an expert, so it gets recomputed a
# lot. (Probably too much; it could be temporarily assocated with each
# expert for the duration of a given vtb. (WWW ??? The fact that this
# isn't organically normalized, but has to be minimaxed to range it,
# is a little disturbing.)
# ==================================================================

expert_level(expert,case) = max(0,min(1.0,feature_distance(expert.spx,case.features)/g_max_features_distance))

ffmt(a)=let s=string(a);string(s[2],s[5],s[8],s[11],s[14]); end # Feature Format

# Choose experts that know about this sort of case.  Also inits the
# algorithm for this case. (%%% The sorting is a bit out of hand
# here. These searches could be done way more efficiently!)

function choose_vtb_experts(case,algorithm_target_certainty)
    global g_experts
    # ==================================================================
    # The "algorithm" is actully one of the experts, chosen to be a
    # close to the target g_algorithm_target_certainty as
    # possible. The theory here is that the algorithm is as good as
    # the best available expert (modulated byt he target certainty).
    # ??? Could also just create a pseudo-expert, and just give it the
    # target certainty.
    # ==================================================================
    sorted_experts = sort(g_experts, by = expert -> abs(expert_level(expert,case)-algorithm_target_certainty),rev=false) 
    the_algorithm=sorted_experts[1]
    # The number of experts reflects the expertise distribution.  If
    # you have a 1.0, then you can stop, otherwise, pick a few. For
    # the moment we just do this in a very ugly way UUU FFF
    sorted_experts = sort(g_experts, by = expert -> expert_level(expert,case),rev=true) 
    top_expertise = expert_level(sorted_experts[1],case)
    # FFF Maybe do this by sampling biased by the top_expertise?
    top_expertise <= 1.0 ? nexperts = 2 : nothing
    top_expertise <= 0.8 ? nexperts = 3 : nothing
    top_expertise <= 0.7 ? nexperts = 4 : nothing
    top_expertise <= 0.6 ? nexperts = 6 : nothing
    return([sorted_experts[i] for i in 1:nexperts],the_algorithm)
end

# g_cstxsdacol caches the treatment sorted aginst the case -- avoids a
# bunch of sorting, and sorting code

global g_cstxsdacol = [] # case_sorted_txs_and_case_overlap_length

function run_a_vtb(name,algorithm_target_certainty)
    global g_cstxsdacol
    case = a_case()
    txsandcol = [[tx, txmatch(case.features,tx.features)] for tx in g_txs]
    g_cstxsdacol = sort(txsandcol, by = x -> x[2], rev=true)
    (experts,the_algorithm) = choose_vtb_experts(case,algorithm_target_certainty) 
    return(vtb(name,case,[expert_consider_a_case(case,expert) for expert in experts],
               the_algorithm,algorithm_target_certainty,expert_consider_a_case(case,the_algorithm),[],[],0.0))
end

# Tx match is just feature_distance + some gaussian noise.
txmatch(casefs,txfs) = round(feature_distance(casefs,txfs)+randn(g_rng,Float64)/5.0; sigdigits = 3)

mutable struct opinion; expert; case; rec; certainty; end

# When an expert considers a case, we start by sorting the treatments
# for this case (feature overlap). Then we sample 5 from the top n
# txs, where n is based on the level of the expert (0.0-1.0). So, for
# example, a complete expert will sample from ~the top 10, etc.

function expert_consider_a_case(case, expert)
    certainty=expert_level(expert,case)
    options=opinion(expert,case,find_treatments(case.features,certainty),certainty)
    return(options)
end

function find_treatments(case_features,certainty)
    #@printf("In FIND_TREATMENTS:\n   certainty=%s\n",certainty)
    selection_range = min(g_ntxs,convert(Int64,round((5/(0.1+(certainty/2.0))))))
    #@printf("   selection_range:%s\n",selection_range)
    weights = Weights([pdf(Poisson(selection_range/5.0),i) for i in 1:selection_range])
    #@printf("   weights:%s\n",weights)
    indexes = [sample(g_txindeces,weights) for i in 1:selection_range]
    #@printf("   indexes:%s\n",indexes)
    samples = g_cstxsdacol[indexes]
    #@printf("   samples:%s\n",samples)
    us = unique(samples)
    return([us[i] for i in 1:min(5,length(us))])
end

global g_all_vtbs = []

function run_vtb_set(prefix;nvtbs=10,algorithm_target_certainty=rand(1:5)/5,verbose=true)
    global g_all_vtbs
    init_vtbs()
    these_vtbs = []
    for v in 1:nvtbs
        vtb = run_a_vtb(string(prefix,"_",length(g_all_vtbs)),algorithm_target_certainty)
        analyze_vtb(prefix,vtb,verbose)
        push!(g_all_vtbs,vtb)
        push!(these_vtbs,vtb)
    end
    return(these_vtbs)
end

global overall_means_of_inter_opinion_distances = [] # Overall Means of Inter-Opinion Distances
global overall_means_of_algorithm_expert_distances = [] # Overall Means of Algorithm-Expert Distances
global overall_means_of_vtb_summary_scores = [] # Overall mean of Overall Scores

function analyze_vtb_set(prefix,vtbs,algorithm_target_certainty)
    global g_all_vtbs, overall_means_of_inter_opinion_distances, overall_means_of_algorithm_expert_distances, overall_means_of_vtb_summary_scores
    @printf("\n\n-----------------------------------\nCurrent VTB set has %s vtbs @ certainty: %s\n",length(vtbs), algorithm_target_certainty)
    omiods = round(mean(map(vtb -> mean(vtb.inter_opinion_distances), vtbs));sigdigits=3)
    omaeds = round(mean(map(vtb -> mean(vtb.algorithm_expert_distances), vtbs));sigdigits=3)
    omovss = round(mean(map(vtb -> vtb.overall_score, vtbs));sigdigits=3)
    @printf("Overall mean of inter-expert distances = %s\n",omiods)
    @printf("Overall mean of algorithm-expert distances = %s\n",omaeds)
    @printf("Overall mean of overall scores = %s\n",omovss)
    push!(overall_means_of_inter_opinion_distances,omiods)
    push!(overall_means_of_algorithm_expert_distances,omaeds)
    push!(overall_means_of_vtb_summary_scores,omovss)
end

txfmt(rec)=@sprintf("[%s:%s]",rec.name,ffmt(rec.features))

function analyze_vtb(prefix,vtb,verbose)
    verbose ? @printf("\nVTB #%s, case features: %s\n",vtb.name,ffmt(vtb.case.features)) : nothing
    top_score = g_cstxsdacol[1][2]
    verbose ? @printf("Top 5 drugs: %s \n",[string("[",txfmt(g_cstxsdacol[i][1]),"@",string(g_cstxsdacol[i][2]),"]") for i in 1:5]) : nothing
    the_algorithm=vtb.the_algorithm
    algorithm_rec=vtb.algorithm_recs.rec
    algorithm_actual_certainty=expert_level(the_algorithm,vtb.case)
    verbose ? @printf("  Algorithm:%s(@%s(%s)):%s\n",ffmt(the_algorithm.spx),algorithm_actual_certainty,
                      vtb.algorithm_target_certainty,[txfmt(rec[1]) for rec in algorithm_rec]) : nothing
    axdists = []
    for opinion in vtb.opinions
        expert = opinion.expert
        txs = [txfmt(rec[1]) for rec in opinion.rec]
        axdist = vtb_tx_distance(algorithm_rec,opinion.rec)
        verbose ? @printf("  Expert #%s:%s:%s(%s)=>%s\n",expert.name,ffmt(expert.spx),txs,opinion.certainty,axdist) : nothing
        push!(axdists,axdist)
    end
    # ==================================================================
    # Scoring inter-opinion v. algorithm-expert scores. At the end of
    # the day we want an overall vtb score that is the mean of the
    # algorithm-expert distances, where each is independently
    # corrected by the "difficulty", which is the mean of the
    # inter-opinion distances. This is all very tricky, in no small
    # part because of the unfortunately terminological choices where
    # "opinion" is an expert's rankings, so you get "inter-opinion"
    # distances, whereas the other distance is called
    # "algorithm-expert" distance. Also, because these are distances,
    # the math seems upside-down, multiplies and divids need to be
    # thought through carefully, and there are places where we add 1
    # to keep things from dividing by zero, etc. Pretty complex for
    # just a few lines of code!
    # ==================================================================
    vtb.inter_opinion_distances=vtb_inter_opinion_dists(vtb.opinions)
    miods = mean(vtb.inter_opinion_distances)
    verbose ? @printf("Inter-expert distances: %s, mean=%s\n", vtb.inter_opinion_distances, round(miods; sigdigits = 2)) : nothing
    vtb.algorithm_expert_distances = axdists
    verbose ? @printf("Algorithm-expert distances: %s, mean=%s\n", axdists, round(mean(axdists); sigdigits = 2)) : nothing
    # 1+ takes care of the case where you have a 0, or very small number, making the overall score essentially infinite (or uncomputable)
    miods < 0.1 ? @printf("\n\n** mean(vtb.inter_opinion_distances) is %s for iods = %s\n\n",miods,vtb.inter_opinion_distances) : nothing
    algorithm_expert_distances_corrected_by_inter_expert_distances = [x/(1+miods) for x in axdists]
    verbose ? @printf("algorithm_expert_distances_corrected_by_inter_expert_distances (+1) = %s\n", algorithm_expert_distances_corrected_by_inter_expert_distances) : nothing
    vtb.overall_score=round(mean(algorithm_expert_distances_corrected_by_inter_expert_distances);sigdigits=3)
    verbose ? @printf("Overall score (mean of just above): %s\n", vtb.overall_score) : nothing
end

function vtb_inter_opinion_dists(opinions)
    txs = map(o->map(t->t[1].name,o.rec),opinions)
    dists = []
    for (index1,r1) in enumerate(txs)
        for (index2,r2) in enumerate(txs[index1+1:length(txs)])
            push!(dists,unequal_hamming(r1,r2))
        end
    end
    return(dists)
end

# Hamming only works on equal-length arrays. We just compare the heads.

function unequal_hamming(r1,r2)
    l1 = length(r1)
    l2 = length(r2)
    if l1 == l2
        return(hamming(r1,r2))
    elseif l1 < l2
        return(hamming(r1,r2[1:l1]))
    else
        return(hamming(r1[1:l2],r2))
    end
end

# Tx_distance for vtbs is more complex than when we knew that the
# target order was always [1,2,3,4,5]. In the VTB world, the lengths
# of the entries can be unequal (although we assume that there are
# always at least 2! (???)).

function vtb_tx_distance(a,b)
    if length(a)==1 || length(b)==1
        @printf("\n\n** VTB_TX_DISTANCE got: \n  a=%s\n  b=%s\n  Returning 25\n\n",a,b)
        return(25)
    end
    # If the first two are the same, or inverted, we ignore everything
    # else and give a 0 or 1.
    if (a[1]==b[1] && a[2]==b[2]); return(0); end
    if (a[1]==b[2] && a[2]==b[1]); return(1); end
    if (a[1]==b[1]); return(2); end
    if (a[1]==b[2] || a[2]==b[1]); return(3); end
    if (a[2]==b[2]); return(4); end
    # Othewise, we just return 2+ the extended hamming distance
    return(4+unequal_hamming(a,b))
end

# This is butt UUUgly but makes it slightly convenient for inputting
# to spreadsheets via c/p w/o having to scp the data back and upload
# it.

function report(label,value,o)
    if typeof(value)<:Dict
        write(o,label); write(o,string(value))
        println(label); println(value)
    else
        println(string(label," = ",value))
        write(o,label)
        s=string("\n",value,"\n")
        write(o,replace(replace(replace(s,r" "=>s"\n"),r"\[|\]|\||\,"=>s""),r"[a-zA-Z]"=>s"")) # UUU !!!
    end
end

function vtb_validation_scan(scan_target_certainties=collect(1:5)/5;verbose=true,nvtbs=10)
    global g_all_vtbs, overall_means_of_inter_opinion_distances, overall_means_of_algorithm_expert_distances, overall_means_of_vtb_summary_scores
    prefix = Dates.format(now(),"yyyymmdd_HHMM_SSssss")
    overall_means_of_inter_opinion_distances = []; overall_means_of_algorithm_expert_distances = []; g_all_vtbs = []; g_all_vtbs = []; overall_means_of_vtb_summary_scores = []
    for algorithm_target_certainty in scan_target_certainties
        @printf("\n\n===============================\nAlgorithm Certainty = %s\n", algorithm_target_certainty)
        analyze_vtb_set(prefix,run_vtb_set(prefix;algorithm_target_certainty=algorithm_target_certainty, verbose=verbose, nvtbs=nvtbs), algorithm_target_certainty)
    end
    open(string("results/",prefix,".xls"), "w") do xls
        report("Parameters",Dict("Total VTBs" => length(g_all_vtbs),
                                 "nvtbs(per scan point)" => nvtbs,
                                 "g_max_features_distance" => g_max_features_distance,
                                 "g_feature_kurtosis" => g_feature_kurtosis,
                                 "g_nexperts" => g_nexperts,
                                 "g_ntxs" => g_ntxs),xls)
        report("Scanned certainties", scan_target_certainties,xls)
        report("Overall means of inter-expert distances",overall_means_of_inter_opinion_distances,xls)
        report("Overall means of algorithm-expert distances",overall_means_of_algorithm_expert_distances,xls)
        report("Overall means of overall scores",overall_means_of_vtb_summary_scores,xls)
    end
end

vtb_validation_scan(;verbose=true,nvtbs=10)
