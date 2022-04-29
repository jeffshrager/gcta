# julia v2a.jl

# ToDo:
# Confirm and document result 1: Linearty of 0.8-1.0
#                             2: no appreciable diff btwn concensus and separate distancing
# Test: Does uniform distance in footrule slow down learning?
#       -> First two positions? (test)
# Antirecmmendations

# ??? Why are the uncorrected summary values LOWER than the corrected
# ones?! The correction should divide by a number that is usually >1!

using Random, Distributions, Printf, StatsBase, Distances, Dates, SciPy
using PyPlot; const plt = PyPlot

# This version is much simpler than valid.jl, based purely on
# distributions, rather than specific feature sets. Cases don't have
# features. We just choose a difficulty for the case, and then create
# a panel of experts, again with a distribution of expertise, based
# upon the case difficulty, and treatments are just indexes, sampled
# from 1:N based on the expertise of each expert (including the
# algorithm).

global g_n_txs = 20
global g_txs = collect(1:g_n_txs)
#@sprintf("g_txs=%s\n",g_txs)

global g_tracedepth = 0
trace!(msg,lvl=0) = lvl >= g_tracedepth ? print(msg) : nothing

function rvtb(algorithm_certainty,distance_function,correction_distance_function)
    trace!(@sprintf("\n\n----- %s ----- Distance function: %s \n", algorithm_certainty, distance_function))
    global g_txs
    # Experts first:
    expert_certainty = min(0.95,1-rand(Beta(1,5))) # 1.0s (or close) screw up the process...and are unrealistic anyway
    trace!(@sprintf("expert_certainty=%s\n",expert_certainty))
    n_experts = 2+convert(Int64,round(6*(1-expert_certainty))) # There's always at least 2 experts
    trace!(@sprintf("n_experts=%s\n",n_experts))
    # COMMENT ME!
    weights=Weights([pdf(Beta(1,2/(1-expert_certainty)),i/g_n_txs) for i in 1:g_n_txs])
    erecs=([sample(g_txs,weights,5;replace=false) for x in 1:n_experts]) # expert recommendations
    trace!(@sprintf("erecs=%s\n",erecs))
    ieds = inter_rec_dists(erecs,correction_distance_function) # inter-expert distances
    trace!(@sprintf("ieds=%s\n",ieds))
    mieds=mean(ieds) # mean of inter-expert distances
    trace!(@sprintf("mieds=%s\n",mieds))
    concensus = find_concensus(erecs)
    trace!(@sprintf("concensus recommendation=%s\n",concensus))
    # Algorithm's turn:
    trace!(@sprintf("algorithm_certainty=%s\n",algorithm_certainty))
    # Alg. rec (Again, -0.01 to avoid divide by 0 @ 1.0 ??? UUU)
    weights=Weights([pdf(Beta(1,2/(1-(algorithm_certainty-.01))),i/g_n_txs) for i in 1:g_n_txs]) 
    arec=sample(g_txs,weights,5;replace=false)
    trace!(@sprintf("arec=%s\n",arec))
    uaedists=[distance_function(arec,erec) for erec in erecs] # Uncorrected algorithm-expert distances
    trace!(@sprintf("uaedists=%s\n",uaedists))
    umaedists=mean(uaedists) # uncrorrected mean of alg-exp distances
    trace!(@sprintf("umaedists=%s\n",umaedists))
    caedists = uaedists/mieds # corrected distances
    trace!(@sprintf("caedists=%s\n",caedists))
    cmaedists = mean(caedists) # mean corrected distances
    trace!(@sprintf("cmaedists=%s\n",cmaedists))
    ucondist=distance_function(arec,concensus)
    trace!(@sprintf("ucondist=%s\n",ucondist))
    ccondist=ucondist/mieds
    trace!(@sprintf("ccondist=%s\n",ccondist))
    # [uncorrected, Corrected, uncorrected_concensus, corrected_concensus]
    return([rn(umaedists;sd=3), rn(cmaedists;sd=3), rn(ucondist;sd=3), rn(ccondist;sd=3)]) 
end
    
function find_concensus(erecs)
    return(erecs[1])
end


function inter_rec_dists(recs,correction_distance_function)
    dists = []
    for (index1,r1) in enumerate(recs)
        for (index2,r2) in enumerate(recs[index1+1:length(recs)])
            push!(dists,correction_distance_function(r1,r2))
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

# Generalized version adds non-maching entries at the end, per https://arxiv.org/pdf/1804.05420.pdf

function generalized_wsf(a,b,wfn=(p->exp(1/p))) # 1/(1+p)
    complete = union(a,b)
    a = cat(a,setdiff(complete,a),dims=1)
    b = cat(b,setdiff(complete,b),dims=1)
    w=[wfn(p) for p in 1:length(a)]
    return (rn(weighted_spearman_footrule_normalized(a, b, w);sd=4))
end

function weighted_spearman_footrule(a, b, w)
    d=0
    p=1
    for e in a
        i=findall(x->x==e,a)[1] # !!! Assumes no replicates, and every element is findable in each!
        j=findall(x->x==e,b)[1]
        d=d+(w[p]*abs(i-j))
        p=p+1
    end
    return(d)
end

function weighted_spearman_footrule_max(a, w)
    n=length(a)
    d=0
    p=1
    for i in a
        d=d+(w[p]*abs(i-(n-i+1)))
        p=p+1
    end
    return(d)
end

function weighted_spearman_footrule_normalized(a, b, w)
    num = weighted_spearman_footrule(a, b, w)
    denom = weighted_spearman_footrule_max(a, w)
    return(num/denom)
end

# global g_distance_functions =
#     [euclidean, sqeuclidean, cityblock, totalvariation, chebyshev, hamming, jaccard, braycurtis, cosine_dist, corr_dist, chisq_dist,
#      gkl_divergence, spannorm_dist, bhattacharyya, hellinger, Distances.meanad, Distances.msd, Distances.rmsd, Distances.nrmsd,
#      rec_distance_12priority_then_5,rec_distance_12priority_then_hamming,nkt,generalized_wsf]

#global g_distance_functions = [weightedtau_py, generalized_wsf, rec_distance_12priority_then_hamming, nkt, euclidean, cityblock, hamming, jaccard]

global g_distance_functions = [generalized_wsf] # [weightedtau_py,generalized_wsf]

function demo_distance_fn(fn)
    testarrays=[[1,2,3,4,5],[5,4,3,2,1],[1,5,2,4,3],[1,1,1,1,1],[1,5,1,5,1],[6,6,6,6,6],[1,6,2,7,3]]
    for a1 in testarrays
        for a2 in testarrays
            @printf("%s(%s,%s)->%s\n",string(fn),a1,a2,fn(a1,a2))
        end
    end
end

####################################################################################
##### YOU ARE DEALING WITH DISTANCE FUNCTIONS, NOT SCORES, SO 0 IS GOOD!!!!!!! #####
####################################################################################

function run(n=1000,experid="experid";learning_sequence=[0.4,0.6,0.8,1.0],correct_p = true)
    @printf("####################################################################################\n")
    @printf("##### YOU ARE DEALING WITH DISTANCE FUNCTIONS, NOT SCORES, SO 0 IS GOOD!!!!!!! #####\n")
    @printf("####################################################################################\n")
    open(string("results/",Dates.format(now(),"yyyymmdd_HHMM_SSssss"),"_",experid,".xls"),"w") do xls
    results_dict = Dict()
        for AE_distance_function in g_distance_functions
            @printf("\n==========================================================\n")
            #demo_distance_fn(AE_distance_function)
            for IE_distance_function in (correct_p ? g_distance_functions : [AE_distance_function]) # If not correcting, just use the same fn
                @printf("\n=== %s/%s ===\n",AE_distance_function,IE_distance_function)
                @printf("[Nb. The reason that the distances are 0.5 is that the stupider the algorithm the closer the distance comes to random.]\n")
                write(xls,@sprintf("\nAE_Distance_Function\tIE_Distance_Function\tLearning_Level\tUNcrtd\tIE_Crtd\tConcensus_UNcrtd\tConcensus_IE_Crtd\n"))
                ccordata = []; ucordata=[] 
                cconcordata = []; uconcordata=[] 
                try
                    for algorithm_certainty in learning_sequence
                        # [uncorrected, Corrected, uncorrected_concensus, corrected_concensus]
                        results = [rvtb(algorithm_certainty,AE_distance_function,IE_distance_function) for i in 1:n]
                        # Separate the results
                        uresults = map(r->r[1],results)
                        cresults = map(r->r[2],results)
                        uconresults = map(r->r[3],results)
                        cconresults = map(r->r[4],results)
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
                        @printf("algcert=%s, mean uncrtd VTB-alg dist=%s (std=%s,n=%s), crtd=%s (%s,%s),uncorr concensus=%s (%s,%s), crtd=%s (%s,%s)\n",
                                algorithm_certainty,
                                muresults,rn(std(uresults);sd=2),length(uresults),
                                mcresults,rn(std(cresults);sd=2),length(cresults),
                                muconresults,rn(std(uconresults);sd=2),length(uconresults),
                                mcconresults,rn(std(cconresults);sd=2),length(cconresults)
                                )
                        write(xls,@sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",AE_distance_function,IE_distance_function,algorithm_certainty,
                                           muresults, mcresults,muconresults, mcconresults))
                        # Save for correlations
                        push!(ucordata,muresults)
                        push!(ccordata,mcresults)
                        push!(uconcordata,muconresults)
                        push!(cconcordata,mcconresults)
                    end
                    # UUU "cor" confusingly both the "cor"relation, as well as meaning "cor"rected 
                    # as: ucor is the uncorrected correlation and ccor is the corrected correlation
                    ucor = cor(ucordata,learning_sequence)
                    ccor = cor(ccordata,learning_sequence)
                    uconcor = cor(uconcordata,learning_sequence)
                    cconcor = cor(cconcordata,learning_sequence)
                    write(xls,@sprintf("%s\t%s\tucor\t%s\tccor\t%s\tuconcor\t%s\tcconcor\t%s\n",AE_distance_function,IE_distance_function,ucor,ccor,uconcor,cconcor))
                    @printf("Correlations between learning and DISTANCEs [nb. negative is good!]: UCor=%s, CCor=%s UConCor=%s, CConCor=%s\n",rn(ucor),rn(ccor),rn(uconcor),rn(cconcor))
                    results_dict[@sprintf("%s/%s",AE_distance_function,IE_distance_function)] = [ucordata,ccordata,uconcordata,cconcordata]
                catch err
                    @printf("******** Failed: %s/%s: %s\n",AE_distance_function,IE_distance_function,string(err))
                end
            end
        end    
        write(xls,"\t\t")
        for elt in learning_sequence
            write(xls,string("\t",elt))
        end
        write(xls,"\n")
        for key in collect(keys(results_dict))
            for normalized in ["normalized","unnormalized"]
                for (label, dx) in zip(["uncorrected","corrected","uncorrected_concensus","corrected_concensus"],1:4)
                    r = results_dict[key][dx]
                    normalized == "normalized" ? normalizer = mean(r) : normalizer = 1.0
                    write(xls,string(normalized,"\t",label,"\t",key))
                    for elt in r
                        write(xls,string("\t",elt/normalizer))
                    end
                    write(xls,"\n")
                end
            end
        end
    end
end

rn(x;sd=2) = round(x;sigdigits=sd)

g_tracedepth = 10
run(1000;learning_sequence=collect(0.1:0.02:1.0),correct_p=false)

