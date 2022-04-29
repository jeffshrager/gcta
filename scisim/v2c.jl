# julia v2c.jl

# Like v2b, c includes contra-recommendations, but the scoring
# algorithm is much better worked out than in previous versions.

# ToDo:
# Confirm and document result 1: Linearty of 0.8-1.0
#                             2: no appreciable diff btwn concensus and separate distancing
# Test: Does uniform distance in footrule slow down learning?
#       -> First two positions? (test)
# Should the Beta paramters for the contrat-recommendations be the formal reverse of the pro?

# ??? Why are the uncorrected summary values LOWER than the corrected
# ones?! The correction should divide by a number that is usually >1!

using Random, Distributions, Printf, StatsBase, Distances, Dates, SciPy, LinearAlgebra
using PyPlot; const plt = PyPlot
include("footrules.jl") # Asher's code

# This version is much simpler than valid.jl, based purely on
# distributions, rather than specific feature sets. Cases don't have
# features. We just choose a difficulty for the case, and then create
# a panel of experts, again with a distribution of expertise, based
# upon the case difficulty, and treatments are just indexes, sampled
# from 1:N based on the expertise of each expert (including the
# algorithm).

global g_n_pros=5
global g_n_cons=2

# This same structure is used for both the pro and con recommendations.
struct distances; erecs; arecs; concensus; uncorrected_individual; corrected_individual; uncorrected_concensus; corrected_concensus; end

global g_n_txs = 20
global g_txs = collect(1:g_n_txs)

global g_tracedepth = 0 #  -1 turns off; Higher numbers trace more (deeper)
trace!(msg,lvl=0) = lvl <= g_tracedepth ? print(msg) : nothing

# This is UUU Ugly: In order to document the function-like parameters
# we need to define string-versions of them. !!!KEEP THESE IN SYNC!!!

borda_score_fn(p) = exp(1/p)
global g_borda_score_fn_string = "exp(1/p)"

global g_max_expert_certainty=0.95

expert_certainty_fn()=1-rand(Beta(1,5))
global g_expert_certainty_fn_string="1-rand(Beta(1,5))"

alphafn(certainty)=1
betafn(certainty)=2/(1-certainty)

global g_alphafn_string="1"
global g_betafn_string="2/(1-certainty)"

function rvtb(algorithm_certainty,distance_function,correction_distance_function)
    global g_n_pros, g_n_cons
    trace!(@sprintf("\n\n===== %s =====\n  (Distance function: %s) \n\n", algorithm_certainty, distance_function))
    # ??? Does the 5 in the beta need to be g_n_pros (which happens to be 5)???
    expert_certainty = min(g_max_expert_certainty,expert_certainty_fn()) # 1.0s (or close) screw up the process...and are unrealistic anyway
    trace!(@sprintf("expert_certainty=%s\n",expert_certainty))
    trace!(@sprintf("algorithm_certainty=%s\n",algorithm_certainty))
    n_experts = 2+convert(Int64,round(6*(1-expert_certainty))) # There's always at least 2 experts
    trace!(@sprintf("n_experts=%s\n",n_experts),2)
    pro_distances=gather_opinions("pro",alphafn,betafn,expert_certainty,algorithm_certainty,
                                  n_experts,distance_function,correction_distance_function, g_n_pros)
    trace!(@sprintf("pro_distances=%s\n",pro_distances),2)
    con_distances=gather_opinions("con",alphafn,betafn,expert_certainty,algorithm_certainty,
                                  n_experts,distance_function,correction_distance_function, g_n_cons)
    trace!(@sprintf("con_distances=%s\n",con_distances),2)
    combined_distances=gather_opinions("combined",alphafn,betafn,expert_certainty,algorithm_certainty,
                                                n_experts,distance_function,correction_distance_function, [g_n_pros,g_n_cons])
    trace!(@sprintf("combined_distances=%s\n",combined_distances),2)
    return(pro_distances,con_distances,combined_distances)
end
    
# This is used for the pros, cons, and combined, just with different
# params. Importantly, if the label is "combined", then n_recs is an
# array: [n_pros,n_cons]. 

function find_options(label,alphafn,betafn,cert,n_reports,n_recs)
    weights=[pdf(Beta(alphafn(cert),betafn(cert)),i/g_n_txs) for i in 1:g_n_txs]
    if label=="con"
        weights=reverse(weights)
    elseif label=="combined"
        weights=cat(weights[1:10],weights[11:20],dims=1)
        n_recs=n_recs[1]+n_recs[2] # Hack to sample the right total number
    end
    trace!(@sprintf("weights=%s\n",weights))
    return [sample(g_txs,Weights(weights),n_recs;replace=false) for x in 1:n_reports]
end

function gather_opinions(label,alphafn,betafn,ecert,acert,n_experts,distance_function,correction_distance_function,n_recs)
    global g_txs, g_n_txs
    trace!(@sprintf("--- %s ---\n",label))
    # Experts first:
    erecs=find_options(label,alphafn,betafn,ecert,n_experts,n_recs)
    trace!(@sprintf("erecs=%s\n",erecs))
    ieds = inter_rec_dists(erecs,correction_distance_function) # inter-expert distances
    trace!(@sprintf("ieds=%s\n",ieds),2)
    mieds=mean(ieds) # mean of inter-expert distances
    trace!(@sprintf("mieds=%s\n",mieds),2)
    concensus = borda_concensus(erecs)
    trace!(@sprintf("concensus %s recommendation=%s\n",label,concensus))
    # Algorithm next:
    arec=find_options(label,alphafn,betafn,acert,1,n_recs)[1] # A bit kludgey UUU
    trace!(@sprintf("%s arec=%s\n",label,arec))
    uaedists=[distance_function(arec,rec) for rec in erecs] # Uncorrected algorithm-expert distances
    trace!(@sprintf("uaedists=%s\n",uaedists),2)
    umaedists=mean(uaedists) # uncrorrected mean of alg-exp distances
    trace!(@sprintf("umaedists=%s\n",umaedists),2)
    caedists = uaedists/mieds # corrected distances
    trace!(@sprintf("caedists=%s\n",caedists),2)
    cmaedists = mean(caedists) # mean corrected distances
    trace!(@sprintf("cmaedists=%s\n",cmaedists),2)
    uconsensusdist=distance_function(arec,concensus)
    trace!(@sprintf("uconsensusdist=%s\n",uconsensusdist),2)
    cconsensusdist=uconsensusdist/mieds
    trace!(@sprintf("cconsensusdist=%s\n",cconsensusdist),2)
    return(distances(erecs, arec, concensus, rn(umaedists;sd=3), rn(cmaedists;sd=3), rn(uconsensusdist;sd=3), rn(cconsensusdist;sd=3)))
end

# There are many ways we could be doing concensus. This is a
# borda-like version, but I use a power function for scoring so that
# the heads get way more weigt. WWW This assumes that the recs are
# all the same length, and uses the length of the first one as the
# template for the result. Thus, any lower-scoring options get dropped
# in thie concensus building page. FFF We may want to return all, and
# trim above.

function borda_concensus(recs)
    len=length(recs[1])
    global g_n_txs
    #g_n_txs=20 # For testing
    # Create a tally sheet
    tallies = [0.0 for i in 1:g_n_txs]
    # Tally up from the end
    for rec in recs
        for p in 1:len
            vote = rec[p]
            score = borda_score_fn(p)
            tallies[vote]+=score
        end
    end
    indexes=findall(p->p>0,tallies)
    return(map(ab->ab[1],sort(map((a,b)->(a,b),indexes,tallies[indexes]),by=ab->ab[2],rev=true))[1:len])
end
            
function inter_rec_dists(recs,correction_distance_function)
    dists = []
    for (index1,r1) in enumerate(recs)
        for (index2,r2) in enumerate(recs[index1+1:length(recs)])
            push!(dists,correction_distance_function(r1,r2,:one))
        end
    end
    return(dists)
end

function rec_distance_12priority_then_hamming(a,b)
    # If the first two are the same, or inverted, we ignore everything
    # else and give a 0 or 1.
    if (a[1]==b[1] && a[2]==b[2]); return(0); end
    if (a[1]==b[2] && a[2]==b[1]); return(1); end
    if (a[1]==b[1]); return(2); end
    if (a[1]==b[2] || a[2]==b[1]); return(3); end
    if (a[2]==b[2]); return(4); end
    # Othewise, we just return 4+ the extended hamming distance
    return(4+hamming(a,b))
end

function rec_distance_12priority_then_5(a,b)
    # If the first two are the same, or inverted, we ignore everything
    # else and give a 0 or 1.
    if (a[1]==b[1] && a[2]==b[2]); return(0); end
    if (a[1]==b[2] && a[2]==b[1]); return(1); end
    if (a[1]==b[1]); return(2); end
    if (a[1]==b[2] || a[2]==b[1]); return(3); end
    if (a[2]==b[2]); return(4); end
    return(5)
end

# Normalized Kt from https://en.wikipedia.org/wiki/Kendall_tau_distance

function nkt(x,y)
    v = 0
    len = length(x)
    for i in 1:len
        for j in (i+1):len
            a = x[i] < x[j] && y[i] > y[j]
            b = x[i] > x[j] && y[i] < y[j]
            (a || b) ? v=v+1 : nothing
        end
    end
    return abs(v)/((len*(len-1)/2.0))
end

function weightedtau_py(x,y)
    len=length(x)
    # Most experiments with different arguments make almost no
    # difference. Except for the weigher function. But in that case a
    # lot of alternate verions of the lambda function that aren't
    # exactly (x->1/(x+1)) end up just returning 1.0. (Probably
    # because they simply hit the ceiling).
    (a,b)=stats.weightedtau(x,y,false) # Normal...same as: (x->1/(x+1)) [Turning off scipy function's lexicographic ranking!]
    #(a,b)=stats.weightedtau(x,y,false,(x->1/(x+2))) # Normal...same as: (x->1/(x+1)) [Turning off scipy function's lexicographic ranking!]
    r = rn((1-((a/2)+0.5));sd=4)
    #println(x,y,r)
    return(r)
end

#wsf_stepfn(p) = exp(1/p) # (1/(1+p))
#global g_wsf_stepfn_string = "exp(1/p)"

global e = exp(1)
wsf_stepfn(p) = (1/(1+p))
global g_wsf_stepfn_string = "UNDEFINED"

# For Asher's generalized_spearman_footrule you need to give it two
# sequences that are each dense in the whole range (in both sense of
# "whole", like "whole numbers", and "entire"), and a "swap cost"
# (which is called position_weights in Asher's code) that contained
# n-1 values.

function wsf(a1,b1,report=:one)
    # These not only have to be the same length, but also have to have
    # the same elements, so if there are missing ones, they get added,
    # but this is asymmetric! 
    a2 = cfr(a1,b1)
    b2 = cfr(b1,a1)
    # Renumber so the vectors are 1-origin and integer dense
    allvals = sort(union(a1,b1),rev=true)
    a3=[findall(p->p==i,allvals)[1] for i in a2]
    b3=[findall(p->p==i,allvals)[1] for i in b2]
    length(a3) != length(b3) && throw("In wfs: a and b must be the same length!")
#    @printf("(%s/%s) -> (%s/%s) -> (%s/%s)\n",a1,b1,a2,b2,a3,b3)
    ews=[wsf_stepfn(p) for p in 1:length(a2)]
    r1 = generalized_spearman_footrule(a3,b3;element_weights=ews)
    r2 = generalized_spearman_footrule(a3,b3;position_weights=ews[1:length(ews)-1])
    r3 = generalized_spearman_footrule(a3,b3;position_weights=ews[2:length(ews)])
    # @printf("  ews=%s\n",ews)
    # @printf("  element_weights=ews :: %s\n",r1)
    # @printf("  position_weights=ews[1:length(ews)-1] :: %s\n",r2)
    # @printf("  position_weights=ews[2:length(ews)] :: %s\n",r3)
    if report == :one
        return(r1) # This is for returning just one score for actual operation (vtbs use this)
    else
        return([r1,r2,r3]) # This is for comparing the weights (demo_distfn uses this mainly)
    end
end

# In cfr (Combine for ranking), a is treated as the target, and b is
# the source. They are assumed to be the same length. (That's checked
# above.) Basically we just add any elementsin b that are not in a, to
# the end of a. 

cfr(a, b)=cat(a,setdiff(b,a),dims=1)

function demo_distfn(distfn)
    rs=[]
    for i in 1:100
        a = sample(1:12, 8, replace = false)
        b = sample(1:12, 8, replace = false)
        r = [a,b,distfn(a,b,:all)]
#        println(r)
        push!(rs,r)
    end
    for entry in sort(rs,by=x->x[3][1])
        println(entry)
    end
end

# This was my old footrule code, but I use Asher's now.
# Generalized version adds non-maching entries at the end, per https://arxiv.org/pdf/1804.05420.pdf
# function generalized_wsf(a,b,wfn=(p->wsf_stepfn(p))) 
#     complete = union(a,b)
#     a = cat(a,setdiff(complete,a),dims=1)
#     b = cat(b,setdiff(complete,b),dims=1)
#     w=[wfn(p) for p in 1:length(a)]
#     return (rn(weighted_spearman_footrule_normalized(a, b, w);sd=4))
# end

# function weighted_spearman_footrule(a, b, w)
#     d=0
#     p=1
#     for e in a
#         i=findall(x->x==e,a)[1] # !!! Assumes no replicates, and every element is findable in each!
#         j=findall(x->x==e,b)[1]
#         d=d+(w[p]*abs(i-j))
#         p=p+1
#     end
#     return(d)
# end

# function weighted_spearman_footrule_max(a, w)
#     n=length(a)
#     d=0
#     p=1
#     for i in a
#         d=d+(w[p]*abs(i-(n-i+1)))
#         p=p+1
#     end
#     return(d)
# end

# Python (from https://arxiv.org/pdf/1804.05420.pdf)

# def weighted_spearman_footrule_max(a,w):
#     n = len(a)
#     d = 0
#     for i, e in enumerate(a,1):
#         # In the original this read: w[e], but that screw python's 0-origin indexing!
#         d += w[e-1] * abs(i - (n - i + 1))
#     return(d)
# weighted_spearman_footrule_max([1, 2, 3, 4, 5, 6, 7],[2.718281828459045, 1.6487212707001282, 1.3956124250860895, 1.2840254166877414, 1.2214027581601699, 1.1813604128656459, 1.1535649948951077])

# function weighted_spearman_footrule_normalized(a, b, w)
#     num = weighted_spearman_footrule(a, b, w)
#     denom = weighted_spearman_footrule_max(a, w)
#     return(num/denom)
# end

# global g_distance_functions =
#     [euclidean, sqeuclidean, cityblock, totalvariation, chebyshev, hamming, jaccard, braycurtis, cosine_dist, corr_dist, chisq_dist,
#      gkl_divergence, spannorm_dist, bhattacharyya, hellinger, Distances.meanad, Distances.msd, Distances.rmsd, Distances.nrmsd,
#      rec_distance_12priority_then_5,rec_distance_12priority_then_hamming,nkt,generalized_wsf]

#global g_distance_functions = [weightedtau_py, generalized_wsf, rec_distance_12priority_then_hamming, nkt, euclidean, cityblock, hamming, jaccard]

global g_distance_functions = [wsf] # [weightedtau_py,wsf]

####################################################################################
##### YOU ARE DEALING WITH DISTANCE FUNCTIONS, NOT SCORES, SO 0 IS GOOD!!!!!!! #####
####################################################################################

function run(n,experid;learning_sequence=[0.4,0.6,0.8,1.0],correct_p = true)
    @printf("####################################################################################\n")
    @printf("##### YOU ARE DEALING WITH DISTANCE FUNCTIONS, NOT SCORES, SO 0 IS GOOD!!!!!!! #####\n")
    @printf("####################################################################################\n")
    date_string=Dates.format(now(),"yyyymmdd_HHMM_SSssss")
    open(string("results/",date_string,"_",experid,".xls"),"w") do xls
        report_params(xls,date_string,experid,n)
        results_dict = Dict()
        for AE_distance_function in g_distance_functions
            @printf("\n==========================================================\n")
            for IE_distance_function in (correct_p ? g_distance_functions : [AE_distance_function]) # If not correcting, just use the same fn
                @printf("\n=== %s/%s ===\n",AE_distance_function,IE_distance_function)
                @printf("[Nb. The reason that the distances are 0.5 is that the stupider the algorithm the closer the distance comes to random.]\n")
                write(xls,@sprintf("\nProCon\tAE_Distance_Function\tIE_Distance_Function\tLearning_Level\tUNcrtd\tIE_Crtd\tConcensus_UNcrtd\tConcensus_IE_Crtd\n"))
                corrdict = Dict("pro"=>Dict("ccordata"=>[],"ucordata"=>[],"cconcordata"=>[],"uconcordata"=>[]),
                                "con"=>Dict("ccordata"=>[],"ucordata"=>[],"cconcordata"=>[],"uconcordata"=>[]),
                                "com"=>Dict("ccordata"=>[],"ucordata"=>[],"cconcordata"=>[],"uconcordata"=>[])
                                )
#                try
                    for algorithm_certainty in learning_sequence
                        distances = [rvtb(algorithm_certainty,AE_distance_function,IE_distance_function) for i in 1:n]
                        # Separate the results
                        pro_distances = map(r->r[1],distances)
                        con_distances = map(r->r[2],distances)
                        com_distances = map(r->r[3],distances)
                        aggregate_and_report("pro",pro_distances,xls,AE_distance_function,IE_distance_function,algorithm_certainty,corrdict)
                        aggregate_and_report("con",con_distances,xls,AE_distance_function,IE_distance_function,algorithm_certainty,corrdict)
                        aggregate_and_report("com",com_distances,xls,AE_distance_function,IE_distance_function,algorithm_certainty,corrdict)
                    end
                    do_correlations("pro",corrdict["pro"],xls,results_dict,learning_sequence,AE_distance_function,IE_distance_function)
                    do_correlations("con",corrdict["con"],xls,results_dict,learning_sequence,AE_distance_function,IE_distance_function)
                    do_correlations("com",corrdict["com"],xls,results_dict,learning_sequence,AE_distance_function,IE_distance_function)
#                catch err
#                    @printf("******** Failed: %s/%s: %s\n",AE_distance_function,IE_distance_function,string(err))
#                end
            end
        end    
        write(xls,"\nCorrelation summaries:\n")
        write(xls,"\t\t")
        write(xls,"\n")
        for key in collect(keys(results_dict))
            for normalized in ["unnormalized","normalized"]
                for (label, dx) in zip(["uncorrected","corrected","uncorrected_concensus","corrected_concensus"],1:4)
                    r = results_dict[key][dx]
                    normalized == "normalized" ? normalizer = mean(r) : normalizer = 1.0
                    write(xls,"\t\t")
                    for elt in learning_sequence; write(xls,string("\t",elt)); end; write(xls,"\n")
                    write(xls,string(normalized,"\t",label,"\t",key))
                    for elt in r; write(xls,string("\t",rn(elt/normalizer,sd=4))); end; write(xls,"\n")
                end
            end
        end
    end
end

function do_correlations(label,dict,xls,results_dict,learning_sequence,AE_distance_function,IE_distance_function)
    # UUU "cor" confusingly both the "cor"relation, as well as meaning "cor"rected 
    # as: ucor is the uncorrected correlation and ccor is the corrected correlation
    ucordata=dict["ucordata"]
    ccordata=dict["ccordata"]
    uconcordata=dict["uconcordata"]
    cconcordata=dict["cconcordata"]
    ucor = cor(ucordata,learning_sequence)
    ccor = cor(ccordata,learning_sequence)
    uconcor = cor(uconcordata,learning_sequence)
    cconcor = cor(cconcordata,learning_sequence)
    write(xls,@sprintf("%s\t%s\t%s\tucor\t%s\tccor\t%s\tuconcor\t%s\tcconcor\t%s\n",
                       label,AE_distance_function,IE_distance_function,ucor,ccor,uconcor,cconcor))
    @printf("%s: Correlations between learning and DISTANCEs [nb. negative is good!]: UCor=%s, CCor=%s UConCor=%s, CConCor=%s\n",
            label,rn(ucor),rn(ccor),rn(uconcor),rn(cconcor))
    results_dict[@sprintf("%s/%s/%s",label,AE_distance_function,IE_distance_function)] = [ucordata,ccordata,uconcordata,cconcordata]
end

function aggregate_and_report(label,distances,xls,AE_distance_function,IE_distance_function,algorithm_certainty,cordict)
    uresults = map(r->r.uncorrected_individual,distances)
    cresults = map(r->r.corrected_individual,distances)
    uconresults = map(r->r.uncorrected_concensus,distances)
    cconresults = map(r->r.corrected_concensus,distances)
    # Delete results that aren't numbers %%% Should use skipmissing, but ... whatever.
    deleteat!(uresults,map(v->(isnan(v)||isinf(v)),uresults)) 
    trace!(@sprintf("Cleaned uncrtd results: %s\n",uresults))
    deleteat!(cresults,map(v->(isnan(v)||isinf(v)),cresults))
    trace!(@sprintf("Cleaned crtd results: %s\n",cresults))
    deleteat!(uconresults,map(v->(isnan(v)||isinf(v)),uconresults))
    trace!(@sprintf("Cleaned un crtd concensus results: %s\n",uconresults))
    deleteat!(cconresults,map(v->(isnan(v)||isinf(v)),cconresults))
    trace!(@sprintf("Cleaned crtd concensus results: %s\n",cconresults))
    # Round to 2 digits
    muresults = rn(mean(uresults);sd=2)
    mcresults = rn(mean(cresults);sd=2)
    muconresults = rn(mean(uconresults);sd=2)
    mcconresults = rn(mean(cconresults);sd=2)
    # Reporting
    @printf("%s: algcert=%s, mean uncrtd VTB-alg dist=%s (std=%s,n=%s), crtd=%s (%s,%s),uncorr concensus=%s (%s,%s), crtd=%s (%s,%s)\n",
            label,algorithm_certainty,
            muresults,rn(std(uresults);sd=2),length(uresults),
            mcresults,rn(std(cresults);sd=2),length(cresults),
            muconresults,rn(std(uconresults);sd=2),length(uconresults),
            mcconresults,rn(std(cconresults);sd=2),length(cconresults)
            )
    write(xls,@sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",label,AE_distance_function,IE_distance_function,algorithm_certainty,
                       muresults, mcresults,muconresults, mcconresults))
    # Save for correlations
    subdict = cordict[label]
    push!(subdict["ucordata"],muresults)
    push!(subdict["ccordata"],mcresults)
    push!(subdict["uconcordata"],muconresults)
    push!(subdict["cconcordata"],mcconresults)
end

function report_params(xls,date_string,experid,n)
    # This would be a LOT simpler if there was something like lisp's
    # symbol-value or if I did something smarter with globals
    global g_n_pros,g_n_cons,g_n_txs,g_borda_score_fn_string,g_wsf_stepfn_string
    write(xls,@sprintf("%s\t%s\n","experid",experid))
    write(xls,@sprintf("%s\t%s\n","time_stamp",date_string))
    write(xls,@sprintf("%s\t%s\n","n",n))
    write(xls,@sprintf("%s\t%s\n","g_n_pros",g_n_pros))
    write(xls,@sprintf("%s\t%s\n","g_n_cons",g_n_cons))
    write(xls,@sprintf("%s\t%s\n","g_n_txs",g_n_txs))
    write(xls,@sprintf("%s\t%s\n","g_borda_score_fn_string",g_borda_score_fn_string))
    write(xls,@sprintf("%s\t%s\n","g_wsf_stepfn_string",g_wsf_stepfn_string))
    write(xls,@sprintf("%s\t%s\n","g_max_expert_certainty", g_max_expert_certainty))
    write(xls,@sprintf("%s\t%s\n","g_expert_certainty_fn_string", g_expert_certainty_fn_string))
    write(xls,@sprintf("%s\t%s\n","g_alphafn_string", g_alphafn_string))
    write(xls,@sprintf("%s\t%s\n","g_betafn_string",g_betafn_string))
end

rn(x;sd=2) = round(x;sigdigits=sd)

g_tracedepth = -1 # -1 turns off; Higher numbers trace more (deeper)
run(1000,"test";learning_sequence=collect(0.8:0.01:0.99),correct_p=false) # Going to 1.0 break the thing!

