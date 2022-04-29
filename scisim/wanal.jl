# julia wanal.jl
# include("wanal.jl")

# Reads from vtb file created by vtb.lisp

using Random, Distributions, Printf, StatsBase, Distances, Dates, SciPy
using CSV, DataFrames, Pingouin, Combinatorics

include("footrules.jl") # Asher's code

# Load up the data and extract case 0 into d0
#d=Pingouin.DataFrame(CSV.File("results/vtb-3855747491.xls"))
#d0 = d[d.p .== 0,:]
# Reshape the results into 5 results from each of the 10 judges
# 5×10 Matrix{Int64}:
#   2   2   2   2   2   2   2   2   2   2
#   3  21   3  28   3   3  21  31  31  31
#   4  41   4  41  21   1  41  41  41  41
#  16  11  41  11   1   4  16  11  16  11
#  11  16  21  42   4  41  11  16  11  16
#dr=reshape(d0[:,:txn],(5,10))

#drt=transpose(dr) # Transpose so that raters are rows
# 10×5 transpose(::Matrix{Int64}) with eltype Int64:
#  2   3   4  16  11
#  2  21  41  11  16
#  2   3   4  41  21
#  2  28  41  11  42
#  2   3  21   1   4
#  2   3   1   4  41
#  2  21  41  16  11
#  2  31  41  11  16
#  2  31  41  16  11
#  2  31  41  11  16

# U will be an array of the union of all selected treatments
#u=reduce(union,dr)
# u=[2, 3, 4, 16, 11, 21, 41, 28, 42, 1, 31]

# Now what we need is an array that is length(u)(=11) high by n
# judges(=10) wide, and is renumbered so that tx 1(=1) is mapped to 1,
# 3->2, 4->3, 16->4, etc., and where whenever a judge does NOT score a
# particular tx, it gets n-selections+1(=6). So, as:

# [[1,1,1,1,1,1,1,1,1,1] # Tx 2 is in position 1 for every judge!
#  [2,6,2,6,2,2,6,6,6,6] # Tx 3 is mostly not included (6), but when it is, it is second
#  [3,6,3,6,5,4,6,6,6,6]
#  [4,5,6,6,6,6,4,5,4,5] # Tx 16 shows up in about half the rankings, mostly in the 4th and 5th positions
#  ...
#  [6,6,6,6,6,6,6,2,2,2]] # That last tx (31) only shows up in the last three judges, but is second in those

function multivtb_concordance(fileid="vtb-3855747491";timestamp,xid="xid",njudges=10,ntoprankings=3,
                              pwsfn=tailharm,metric=ssfr) 
    @printf("\n\n======= multivtb_concordance (%s) [xid=%s] pwsfn=%s, metric=%s =======\n\n",fileid,xid,pwsfn,metric)
    results=[]
    d=Pingouin.DataFrame(CSV.File(string("results/",fileid,".xls")))
    for pt in unique(d.p)
        d0 = d[d.p .== pt,:]
        dr=reshape(d0[:,:txn],(5,njudges)) # Initially there are 5
        # make rows = raters (=rankings) in correct order, top n (3:
        # ignoring cons, and excess pros)) 
        drt=transpose(dr[1:ntoprankings,:])
        mjc=multijudge_concordance(drt;pwsfn=pwsfn,metric=metric)
        con=concensus(drt)
        ace=algorithm_concordance_estimate(con,drt;pwsfn=pwsfn,metric=metric)
        push!(results,[mjc,pt,drt,ace,con])
    end
    open(string("results/",timestamp,"_",xid,"_multicon_",fileid,"_",pwsfn,"_",metric,".xls"), "w") do io
        @printf(io,"PatintNo\tConcordancet\tAlgDistMean[vsJ1]\tAlgDistStD\tConcensus\tJudges\n")
        for e in sort!(results,by=x->x[1],rev=true)
            @printf(io,"%s\t%s\t%s\t%s\t%s",e[2],e[1],e[4][1],e[4][2],e[5])
            for j in 1:njudges
                @printf(io,"\t%s",e[3][j,:])
            end
            @printf(io,"\n")
        end
    end
end

# Creates a concensus and then finds the distance between that and all
# the rest and take the mean and deviation.

function algorithm_concordance_estimate(con,drt;pwsfn=pwsfn,metric=metric)
    dsts=[]
    for i in 1:size(drt)[1]
        c = open_concordance(con,drt[i,:];pwsfn=pwsfn,metric=metric)
        push!(dsts,c)
    end
    return([mean(dsts),std(dsts)])
end

# Finds the concensus of a bunch of judgements. We use a sligtly
# tricky was of scoring for position so that thing that are in the
# first position a lot end up first.

function concensus(drt)
    d=Dict()
    rl = size(drt[1,:])[1]
    for i in 1:size(drt)[1]
        r = drt[i,:]
        for j in 1:rl
            tx = r[j]
            if haskey(d,tx)
                d[tx]=d[tx]+1+(rl-j) # Decrement by position 
            else
                d[tx]=1+(rl-j)
            end
        end
    end
    return(map(x->x[1],sort(collect(d),by=x->x[2],rev=true))[1:rl])
end    

# This requires an array where each row is the rankings of n items
# from one judge. All judges must report the same n (number of items),
# but the items to not have to be the same in each ranking.

function multijudge_concordance(drt;pwsfn=tailharm,metric=ssfr)
    njudges=size(drt)[1]
    sum=0; n=0
    for j1 in 1:njudges
        for j2 in (j1+1):njudges
            r1=drt[j1,:];r2=drt[j2,:]
            gsf=open_concordance(r1,r2;pwsfn=pwsfn,metric=metric)
            sum=sum+gsf
            n=n+1
        end
    end
    return(sum/n)
end

# This sort of concordance takes into account both rank order
# (basically the same as the footrule) but allows for the ranked lists
# to have different elements. Elements that are not listed in one or
# the other are bascially added to the tail end.

# for [1,2,3]x[1,12,13] => [[1,2,3,12,13], [1,2,3,13,12], [1,12,13,2,3] , [1,12,13,3,2]]
# for [1,2,3]x[1,2,13] => [[1,2,3,13], [1,2,13,3]]
# for [1,2,3]x[11,12,13] => [[1,2,3,11,12,13],[1,2,3,11,13,12]+4 more
#                             [11,12,13,1,2,3] .... 5 more]                

function open_concordance(a,b;pwsfn=tailharm,metric=ssfr)
    u=union(a,b)
    abplus=append!([append!(copy(a),more) for more in collect(permutations(setdiff(u,a)))],
                   [append!(copy(b),more) for more in collect(permutations(setdiff(u,b)))])
    lenabplus = length(abplus)
    pws=pwsfn(a,b,abplus)
    # @printf("a=%s, b=%s\n",a,b)
    # @printf("abplus=%s\n",abplus)
    # @printf("pws=%s\n",pws)
    sum=0;n=0
    for i in 1:lenabplus
        for j in (i+1):lenabplus
            ar = abplus[i]; br=abplus[j]
            r=metric(ar,br,pws)
            sum=sum+r; n=n+1
        end
    end
    return(sum/n)
end

function mvtbcdistest(l;timestamp,pwsfn=tailharm,metric=:symftrule,xid="xid")
    @printf("\n\n======= mvtbcdistest (%s) [xid=%s]  pwsfn=%s, metric=%s  =======\n\n",l,xid,pwsfn,metric)
    d=Dict()
    open(string("results/",timestamp,"_",xid,"_mvtbcdistest_",l,"_",pwsfn,"_",metric,".xls"), "w") do io
        for p in 1:10000
            a=nrep(l,2*l);b=nrep(l,2*l)
            c=open_concordance(a,b;pwsfn=pwsfn,metric=metric)
            ip=idperm(a,b)
            d[ip]=c
            @printf(io,"%s\t%s\t%s\t\"%s\"\n",a,b,c,ip)
        end
    end
    open(string("results/",timestamp,"_",xid,"_idperms_",l,"_",pwsfn,"_",metric,".xls"), "w") do io
        for elt in sort(collect(d),by=(x->x[2]))
            @printf(io,"%s\t%s\n",elt[1],elt[2])
            println(elt)
        end
    end
end

function idperm(a,b)
    r=string("x",findornull(a,1,b))
    for i in 2:length(a)
        r=string(r,findornull(a,i,b))
    end
    return(r)
end

function findornull(a,p,b)
    r = findall((x)->x==a[p],b)
    if 0==length(r);return("-");end
    return(r[1])
end

function nrep(len,max)
    r = []
    n = [x for x in 1:max]
    for i in 1:len
        push!(r,rand(setdiff(n,r)))
    end
    return(convert(Vector{Int64}, r))
end

# Whereas Asher's footrule treats the elements in the array as
# indexes, this treats them merely as symbols. Indeed, this will work
# just as well with names, like [:a, :b, :c] or ["a", "b", "c"] as
# with numbers. Things that don't match between ranks are scored as
# the length of the array+1. My position weights are more like Asher's
# element weights. They are multipied in when something has to be
# moved out of a position in the a list, into one in the b list. 

# Note that if the pws are all 1 (or all the same), the symmetric part
# is reduntant. 

function ssfr(a,b,pws) # symmetric_symbolic_footrule
    return((symbolic_footrule(a,b,pws)+symbolic_footrule(b,a,pws))/2.0)
end

function symbolic_footrule(a,b,pws)
    sum=0
    # The miss isn't relevant in the presencee of the concordance
    # because that fills out the tail with the permultaions of missing
    # elt, so there are no missing elts in that case, but the missing
    # score is relevant if this is called with missing values.
    missx=1+length(a)
    # For each elelment in a find its index in b and add up the distance differnce
    j=1
    for elt in a
        pos=findall(x->x==elt,b)
        if pos==[]
            sum=sum+missx*pws[j]
        else
            sum=sum+abs(pos[1]-j)*pws[j]
        end
        j=j+1
    end
    return(sum)
end

# This version of the footrule simply counts the number of times
# things are in the wrwong order. When something does't occur at all,
# you get two points!

function ltgt(a,b,pws)
    return((oneway_ltgt(a,b,pws)+oneway_ltgt(b,a,pws))/2.0)
end

function oneway_ltgt(a,b,pws)
    # Need to * pws into missing
    sum=0
    j=0
    for elt in a
        posb=findall(x->x==elt,b)
        posa=findall(x->x==elt,a)
        if posb==[]
            sum=sum+2
        elseif posb[1]<posa[1] # wrong order!
            sum=sum+pws[j]
        end
        j=j+1
    end
    return(sum)
end

# PWS Algorithms
function all1pws(a,b,abplus)
    return([1.0 for i in 1:length(abplus[1])])
end

function tailharmpws(a,b,abplus) 
    return(append!([1.0/(2^i) for i in 1:length(a)],[1/(2^(2+length(a))) for i in 1:(1+length(abplus[1])-length(a))]))
end

function tailharmpwswithleadning1(a,b,abplus) 
    return(append!([1.0],append!([1.0/(2^i) for i in 1:(length(a)-1)],[1/(2^(2+length(a))) for i in 1:(1+length(abplus[1])-(length(a)+1))])))
end

function randpws(a,b,abplus)
    return([rand(1)[1] for i in 1:length(abplus[1])])
end

# Exploring how the symmetric scales up with the length of the record
function explore1()
    timestamp=Dates.format(now(),"yyyymmdd_HHMM_SSssss")
    for pwsfn in [tailharmpws, all1pws, randpws]
        for metric in [ssfr ltgt]
            for l in 3:5
                mvtbcdistest(l;timestamp=timestamp,pwsfn=pwsfn,metric=metric,xid=string(l,"ex1a"))
            end
            multivtb_concordance("vtb-3855747491";timestamp=timestamp,pwsfn=pwsfn,metric=metric,xid="ex1vtbs")
        end
    end
end

# Exploring different versions of the position weights
function explore2()
    timestamp=Dates.format(now(),"yyyymmdd_HHMM_SSssss")
    for pwsfn in [all1pws, tailharmpws, randpws]
        for l in 3:5
            mvtbcdistest(l;timestamp=timestamp,pwsfn=pwsfn,metric=ssfr,xid=string(l,"ex2"))
        end
    end
end

# Focus on tailharm to make sure that the algorithm isn't wrong
function explore3()
    timestamp=Dates.format(now(),"yyyymmdd_HHMM_SSssss")
    for pwsfn in [tailharmpws, tailharmpwswithleadning1]
        for l in 3:5
            mvtbcdistest(l;timestamp=timestamp,pwsfn=pwsfn,metric=ssfr,xid=string(l,"ex3"))
        end
    end
end


# Focus on tailharm to make sure that the algorithm isn't wrong
function explore3()
    timestamp=Dates.format(now(),"yyyymmdd_HHMM_SSssss")
    for pwsfn in [tailharmpws, tailharmpwswithleadning1]
        for l in 3:5
            mvtbcdistest(l;timestamp=timestamp,pwsfn=pwsfn,metric=ssfr,xid=string(l,"ex3"))
        end
    end
end

# Using tailharmpws.ssfr, compute hold-1-out algorithm contrasts for each judgement
function explore4()
    timestamp=Dates.format(now(),"yyyymmdd_HHMM_SSssss")
    multivtb_concordance("vtb-3855747491";timestamp=timestamp,pwsfn=tailharmpws,metric=ssfr,xid="ex4")
end

#explore1()
#explore2()
#explore3()
explore4()
