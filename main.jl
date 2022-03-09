
## HELPER functions ####################
function combinerraw!(df,name)
    # s=split(name,"_")
    cname = name # string(s[1:end-1]...)
    vname = string(name,"_Vardict_Raw")
    lname = string(name,"_Lofreq_Raw")
    mname = string(name,"_Mutect_Raw")
    df[:,cname] = Any[ismissing(i) ? missing : i for i in df[:,vname]]
    if occursin("IsComplex",cname)
        for i in 1:nrow(df)
            if df[i,mname] + df[i,vname] + df[i,lname] > 0 #(! ismissing(df[i,vname]) && df[i,vname]) || (! ismissing(df[i,mname]) && df[i,mname]) || (! ismissing(df[i,lname]) && df[i,lname])
                df[i,cname] = true
            # else
            #     df[i,cname] = false
            end
        end
    elseif occursin("ComplexKey",cname)
        for i in 1:nrow(df)
            if df[i,mname] != ""
                df[i,cname] = df[i,mname]
            elseif df[i,vname] != ""
                df[i,cname] = df[i,vname]
            elseif df[i,lname] != "" #(! ismissing(df[i,vname]) && df[i,vname]) || (! ismissing(df[i,mname]) && df[i,mname]) || (! ismissing(df[i,lname]) && df[i,lname])
                df[i,cname] = df[i,lname]
            else
                df[i,cname] = ""
            end
        end
    else
        for i in 1:nrow(df)
            if ismissing(df[i,cname]) && ! df[i,"IsComplex"]
                if ! ismissing(df[i,lname])
                    df[i,cname] = df[i,lname]
                elseif ! ismissing(df[i,mname])
                    df[i,cname] = df[i,mname]
                end
            end
        end
    end
    select!(df,Not([lname,vname,mname]))
    return nothing
    
end

function combiner!(df,name)
    #_combine
    cname = string(split(name,"_")[1:end-1]...)
    vname = string(name,"_Vardict")
    lname = string(name,"_Lofreq")
    mname = string(name,"_Mutect")
    df[:,cname] = Any[ismissing(i) ? missing : i for i in df[:,vname]]
    if cname == "IsComplex"
        for i in 1:nrow(df)
            if df[i,vname] + df[i,mname] + df[i,lname] > 0 #(! ismissing(df[i,vname]) && df[i,vname]) || (! ismissing(df[i,mname]) && df[i,mname]) || (! ismissing(df[i,lname]) && df[i,lname])
                df[i,cname] = true
            # elseif ismissing(df[i,cname])
            #     df[i,cname] = false
            end
        end
    else
        for i in 1:nrow(df)
            if ismissing(df[i,cname]) && ! df[i,"IsComplex"]
                if ! ismissing(df[i,lname])
                    df[i,cname] = df[i,lname]
                elseif ! ismissing(df[i,mname])
                    df[i,cname] = df[i,mname]
                end
            end
        end
    end
    select!(df,Not([lname,vname,mname]))
    return nothing
    
end

function sntag!(v)
    for i in 1:length(v)
        v[i] = split(v[i],"__")[1]
    end
end

fpfilters = [
    "RLD25",
    "MVF0",
    "SB1",
    "NRC",
    "PB10",
    "IRC",
    "MMQSD50",
    "DETP20",
    "MVC4",
    "MMQS100",
    "MQD30"
    ]
function test(df, col)
    spl = split.(df[:,col],";")
    fpfil = [ "" for i in 1:length(spl) ]
    caller = [ "" for i in 1:length(spl) ]
    FPpass = [ true for i in 1:length(spl) ]
    for sampl in 1:length(spl)
        for i in 1:length(spl[sampl])

            # println( spl[sampl][i] )

            if spl[sampl][i] in fpfilters
                if spl[sampl][i] != "MVF0" && spl[sampl][i] != "IRC"
                    fpfil[sampl] *= ";" * spl[sampl][i]
                    FPpass[sampl] = false
                end
            else
                # println(caller[sampl],spl[sampl][i])
                caller[sampl] *= ";" * spl[sampl][i]
            end
        end
        if fpfil[sampl] == "" 
            fpfil[sampl] *= "PASS"
        else
            fpfil[sampl] = fpfil[sampl][2:end]
        end
        if caller[sampl] == "" 
            caller[sampl] *= "PASS"
        else
            caller[sampl] = caller[sampl][2:end]
        end
    end

    return fpfil,caller,FPpass
end

function altchg(alt,ref)
    if (alt == "A" && ref == "G") || (alt == "T" && ref == "C")
        return "AGTC"
    elseif (alt == "A" && ref == "T") || (alt == "T" && ref == "A")
        return "ATTA"
    elseif (alt == "C" && ref == "G") || (alt == "G" && ref == "C")
        return "CGGC"
    elseif (alt == "C" && ref == "T") || (alt == "G" && ref == "A")
        return "CTGA"
    elseif (alt == "A" && ref == "C") || (alt == "T" && ref == "G")
        return "ACGT"
    elseif (alt == "C" && ref == "A") || (alt == "G" && ref == "T")
        return "CAGT"
    elseif (alt=="I" || ref=="I")
        return "I"
    else
        error("No good base")
    end
end

function textonehot(df,col,delim=";")
	filtermut = [split(i,delim) for i in df[:,col]]
	filtermutect = unique(reduce(vcat,filtermut))
	filtermutectdata = DataFrame(zeros(Float64,nrow(df),length(filtermutect) ),string.(col,"_",filtermutect) )
	for i in 1:length(filtermut)
		for j in 1:length(filtermut[i]) # j is column to be set to 1
			id=findfirst(x->x == filtermut[i][j],filtermutect)
			filtermutectdata[i,id] = 1.0
		end
	end
	oldnames=names(filtermutectdata)
	nams = replace.(names(filtermutectdata),"."=>"_")
	diffs = oldnames .== nams
	for i in 1:length(nams)
		if ! diffs[i]
			DataFrames.rename!(filtermutectdata,oldnames[i]=>nams[i])
		end
	end
	return filtermutectdata
end

################################################################ SQO #################################

function query(P)

    # return 0
    # rename!(P,:Column2 => :CHROM)
    # rename!(P,:Column3 => :POS)
    # rename!(P,:Column4 => :REF)
    # rename!(P,:Column5 => :ALT)
    # rename!(P,:Column7 => :gt_AD)
    # altref = split.(P.gt_AD,",")
    # altref = P.gt_AD
    P.altcount = copy(P.gt_AD_alt_vardict)
    P.refcount = copy(P.gt_AD_ref_vardict)
    filter!(x->x.altcount < x.refcount, P)
    P.sample_tag = ["sample" for i in P.altcount]
    filter!(x->x.altcount > 1, P)
    P.emp_vaf = P.altcount ./ (P.refcount .+ P.altcount)
    P.passed = trues(nrow(P))
    for i in 1:nrow(P)
        if length(P.REF[i]) == 1 == length(P.ALT[i])
            P.passed[i] = false
        elseif length(P.REF[i]) == 1 && length(P.ALT[i]) < 3
            P.passed[i] = false
        elseif length(P.REF[i]) < 3 && length(P.ALT[i]) == 1
            P.passed[i] = false
        end
    end
    filter!(x->x.passed,P)
    select!(P,Not(:passed))
    # filter!(x-> ! (1 == length(x.REF) && 1 == length(x.ALT)) , P)
    P.Column1 = [i for i in 1:length(P[!,1])]
    P.POS .-= 1
    P.REF = "N" .* P.REF
    P.ALT = "N" .* P.ALT
    return P
end

function subject(M)
    # delm=","
    # if name[end-2:end] == "tsv"
    #     delm = "\t"
    # end
    # M = CSV.read(name,DataFrame,header=1,delim=delm; stringtype=String)
    if sum(occursin.("ofreq",names(M))) > 0
        rename!(M, [:gt_AD_alt_lofreq=>:gt_AD_alt,:gt_AD_ref_lofreq=>:gt_AD_ref])
    elseif sum(occursin.("utect",names(M))) > 0
        rename!(M, [:gt_AD_alt_mutect=>:gt_AD_alt,:gt_AD_ref_mutect=>:gt_AD_ref])
    elseif sum(occursin.("ardict",names(M))) > 0
        rename!(M, [:gt_AD_alt_vardict=>:gt_AD_alt,:gt_AD_ref_vardict=>:gt_AD_ref])
    end
    M.altcount = copy(M.gt_AD_alt)
    M.refcount = copy(M.gt_AD_ref)
    M.sample_tag = ["sample" for i in M.altcount]
    M.emp_vaf = M.altcount ./ (M.refcount .+ M.altcount)
    # filter!(x->x.emp_vaf > 0.00001, M)
    # filter!(x->x.altcount < x.refcount, M)
    M.Column1 = [i for i in 1:length(M[!,1])]
    # M.sample_key = string.(M[:,"SAMPLE"]," ",M[:,"CHROM"]," ",M[:,"POS"]," ",M[:,"REF"],">",M[:,"ALT"])
    # M.key = string.(M[:,"CHROM"]," ",M[:,"POS"]," ",M[:,"REF"],">",M[:,"ALT"])
    return M
end

function part1(df,x)
    finalgroups2 = []
    perm = sortperm(x.POS) # sort by position
    x = x[perm,:] # reorder by POS
    groups = [Vector{Int32}() for i in 1:length(x.POS)]

    for row in 1:length(x.POS) # go through each variant
        temp = Vector{Int}()
        rww = row
        while true
            rww+=1
            if rww <= length(x.POS) && -1 < x.POS[rww] - x.POS[row] < 20 # see if neighbor varaints are <= 15 bp away
                push!(temp, x.Column1[rww])
            else
                break
            end
        end
        rww = row
        while true
            rww -= 1
            if rww > 0 && -1 < x.POS[row] - x.POS[rww] < 20 # see if neighbor varaints are <= 15 bp away
                push!(temp, x.Column1[rww])
            else
                break
            end
        end

        if 21 > length(temp) > 0 # if group of at least
            push!(temp, x.Column1[row])
            groups[row] = temp
        end
    end
    groups = [m for m in groups if ! isempty(m)]
    if length(groups) > 0 # see if there are any potential complex variants
        for re in 1:length(groups)
            groups[re] = sort(groups[re]) # sort each potential group by index value
        end
        if length(groups) > 1
            grps = [Vector{Int32}() for i in 1:length(groups)] # Vector{Vector{Int32}}()
            cur = groups[1]
            loc = 2
            grps[1] = groups[1]
            for g in 2:length(groups) # check if there are duplicated complex
                if g == length(groups)
                    grps[loc] = groups[g]
                elseif groups[g] != cur
                    grps[loc] = groups[g]
                    loc +=1
                    cur = groups[g]
                end
            end
            grps = [m for m in grps if ! isempty(m)]
        end
        groups = [Vector{Int32}() for i in 1:length(grps)] # Vector{Vector{Int32}}()
        tru = fill(true,length(grps))
        for g in 1:length(grps) # check if any groups are subsets of any others
            for gg in 1:length(grps)
                if tru[gg] && g != gg && issubset(grps[g], grps[gg])
                    tru[g] = false
                    @goto subsetend
                end
            end
            groups[g] = grps[g]
            @label subsetend
        end
        groups = [m for m in groups if ! isempty(m)]
        tempg = Vector{Vector{Int}}()
        # emp = df[:,:emp_vaf]
        # poscount = df[:,:altcount]
        # totalcount = df[:,:refcount]
        for g in groups
            mxx = 6 #length(g) > 6 ? 6 : length(g)
            temp = vcat([collect(Combinatorics.combinations(g,i)) for i in 1:mxx]...) #Vector{Int32}() ## gg is each starting point
            if length(temp) > 0
                append!(tempg,temp)
            end

        end
        if length(tempg) > 0
            for i in tempg
                sort!(i)
            end
            if length(tempg) > 1
                grps = [Vector{Int32}() for i in 1:length(tempg)] # Vector{Vector{Int32}}()
                cur = tempg[1]
                loc = 1
                for g in 2:length(tempg) # cut repeat small groups

                    if g == length(tempg)
                        grps[loc] = tempg[g]
                    elseif tempg[g] != cur
                        grps[loc] = cur
                        loc += 1
                        cur = tempg[g]
                    end
                end
                grps = [m for m in grps if ! isempty(m)]
            end

            groups = grps
            grps = [] #  Vector{Vector{Int32}}()
            alt = df[:,:ALT]
            pos = df[:,:POS]
            reff = df[:,:REF]
            for g in groups # group

                strt = pos[g[1]]
                # if pos[g[end]] + maximum(length.( reff[g] ) ) - strt < 1 || pos[g[end]] + maximum(length.( alt[g] ) ) - strt < 1
                #     @goto badspot
                # end
                complexalt = fill('-',pos[g[end]] + maximum(length.( alt[g] ) ) - strt)
                complexref = fill('-',pos[g[end]] + maximum(length.( reff[g] ) ) - strt)
                complex = fill('-',pos[g[end]] + maximum(length.( alt[g] ) ) + 30 - strt)
                #
                for gg in 1:length(g)  # element

                    for s in 1:length(alt[g[gg]]) # base
                        if complexalt[pos[g[gg]]-strt+s] == '-' || complexalt[pos[g[gg]]-strt+s] == alt[g[gg]][s]
                            complexalt[pos[g[gg]]-strt+s] = alt[g[gg]][s]
                        else
                            @goto badspot
                        end
                    end
                end
                while true
                    if complexalt[end] == '-'
                        pop!(complexalt)
                    else
                        break
                    end
                end

                for gg in 1:length(g)  # element

                    for s in 1:length(reff[g[gg]]) # base
                        if complexref[pos[g[gg]]-strt+s] == '-' || complexref[pos[g[gg]]-strt+s] == reff[g[gg]][s]
                            complexref[pos[g[gg]]-strt+s] = reff[g[gg]][s]
                        else
                            @goto badspot
                        end
                    end
                end
                while true
                    if complexref[end] == '-'
                        pop!(complexref)
                    else
                        break
                    end
                end
                complex = complexalt
                push!(grps,(g,complexalt,length(complexalt),complexref, length(complexref),complex, length(complex)  ))

                @label badspot
            end
            if length(grps) > 0
                push!(finalgroups2,((x.sample_tag[1],x.CHROM[1]),grps))
            end
        end


    end
    return finalgroups2
end

function part2(df,finalgroup,pp)

    temp=[]
    tempfinal = []

    try
        pp[( sample_tag=finalgroup[1][1], CHROM = finalgroup[1][2])]
    catch
        # println("\n Not found \n",finalgroup[1][1])
        # println(finalgroup[1][2])
        return []
    end
    p = pp[( sample_tag=finalgroup[1][1], CHROM = finalgroup[1][2])]

    poss = df[:,:POS]
    emp = df[:,:emp_vaf]
    reff = df[:,:REF]
    altt = df[:,:ALT]
    for j in 1:length(p[!,1])

        for k in 1:length(finalgroup[2]) # each group
           
            # println(length(p.ALT[j]))
            if  0 <= poss[finalgroup[2][k][1][1]] - p.POS[j] && poss[finalgroup[2][k][1][end]] - p.POS[j] <= max(length(p.REF[j]),length(p.ALT[j]))  # check if pindel is larger than complex and within 15 bp of start of complex
                # push!(temp, (p.Column1[j], p.REF[j], p.ALT[j], finalgroup[2][k][6], finalgroup[2][k][4], finalgroup[2][k][2], finalgroup[2][k][1], coninf[j], p.emp_vaf[j], emp[finalgroup[2][k][1]], p.CHROM[j], p.sample_tag[j], p.POS[j], poss[finalgroup[2][k][1]], reff[finalgroup[2][k][1]], altt[finalgroup[2][k][1]], p.altcount[j], string(finalgroup[2][k][6]...) ) )
                push!(temp, (p.Column1[j], p.REF[j], p.ALT[j], finalgroup[2][k][6], finalgroup[2][k][4], finalgroup[2][k][2], finalgroup[2][k][1], [], p.emp_vaf[j], emp[finalgroup[2][k][1]], p.CHROM[j], p.sample_tag[j], p.POS[j], poss[finalgroup[2][k][1]], reff[finalgroup[2][k][1]], altt[finalgroup[2][k][1]], p.altcount[j], string(finalgroup[2][k][6]...) ) )
                #               1            2           3           4                   5                        6                 7                 8            9           10                        11           12              13          14                           15                          16                          17             18                           19
            end
        end
    end
    for as in temp
        
        f = 0
        doublefail = 0
        for i in 1:length(as[6]) # go through and make sure complex is completly within pindel and in order; in ALT
            if as[6][i] != '-'
                a=findfirst(x->x == as[6][i], as[3][f+1:end])

                if isnothing(a)
                    a=findfirst(x->x == 'N', as[3][f+1:end])

                    if isnothing(a)

                        doublefail +=1
                        @goto skip1
                    else
                        f += a
                    end
                else
                    f += a
                end
            end
        end
        @label skip1
        f = 0
        for i in 1:length(as[6]) # go through and make sure complex is completly within pindel and in order; in ALT
            if i == 1
                f+=1
            elseif as[6][i] != '-'
                a=findfirst(x->x == as[6][i], as[3][f+1:end])

                if isnothing(a)
                    a=findfirst(x->x == 'N', as[3][f+1:end])

                    if isnothing(a)

                        doublefail +=1
                        @goto skip2
                    else
                        f += a
                    end
                else
                    f += a
                end
            end
        end
        @label skip2
        if doublefail == 2
            @goto notcomplex
        end
        f = 0
        doublefail = 0
        for i in 1:length(as[5]) # go through and make sure complex is completly within pindel and in order; in REF
            if as[5][i] != '-'
                a=findfirst(x->x == as[5][i], as[2][f+1:end])
                if isnothing(a)
                    a=findfirst(x->x == 'N', as[2][f+1:end])
                    if isnothing(a)

                        doublefail +=1
                        @goto skip3
                    else
                        f += a
                    end
                else
                    f += a
                end
            end
        end
        @label skip3
        f = 0
        for i in 1:length(as[5]) # go through and make sure complex is completly within pindel and in order; in REF
            if i == 1
                f+=1
            elseif as[5][i] != '-'
                a=findfirst(x->x == as[5][i], as[2][f+1:end])
                if isnothing(a)
                    a=findfirst(x->x == 'N', as[2][f+1:end])
                    if isnothing(a)

                        doublefail +=1
                        @goto skip4
                    else
                        f += a
                    end
                else
                    f += a
                end
            end
        end
        @label skip4
        if doublefail == 2
            @goto notcomplex
        end
        push!(tempfinal, as)
        @label notcomplex
    end
    if length(tempfinal) == 0
        return []
    end
    compl = [i[4] for i in tempfinal]
    perm = sortperm(compl)
    return tempfinal[perm]
end

function groupsfun(df,P)
    s=groupby(df, [:sample_tag,:CHROM]) # group by sample_tag and chromosome
  
    finalgroups = [[] for i in 1:length(s)]
    finalgroups[1] = part1(df,s[1])
    for il in 2:length(s) # go through each group
        finalgroups[il] = part1(df,s[il])
    end
    finalgroups = vcat([m for m in finalgroups if ! isempty(m)]...)
    
    final = [[] for i in 1:length(finalgroups)]
    if length(final) == 0
        # println("fail 2")
        return nothing
    else
        for ii in finalgroups
            if length(ii) < 1
                error("ii finalgrousp")
            end
        end
        # println("length part 1 ",length(final))
    end
    pp=groupby(P, [:sample_tag,:CHROM])
    final[1] = part2(df,finalgroups[1],pp)
    # println(finalgroups[1][2][1])
    for i in 2:length(finalgroups)
        final[i] = part2(df,finalgroups[i],pp)
    end
    vv = vcat([m for m in final if ! isempty(m)]...)
    
    if length(vv) == 0
        return nothing
    end
    d=DataFrame(vv)
    # println("part 2 done")
    d.literal = d[:,18]
    d.truncated = d[:,18]
    xx = d
    occurslit = zeros(Bool,length(xx[!,18]))
    occurstrunc = zeros(Bool,length(xx[!,18]))
    for i in 1:length(xx[!,18])
        xx[i,:literal] = replace(xx[i,:literal],"-"=> ".")
        if occursin(Regex(xx[i,:literal]),xx[i,3])
            occurslit[i] = true
        end
        gg = xx[i,:truncated]
        newgg = ""
        for g in 1:length(gg)
            if ! (gg[g] in ['-','.','*'])
                newgg *= gg[g]
            elseif gg[g] == '-' && ! (gg[g-1] in ['-','.','*'])
                newgg *= '.'
            elseif gg[g] == '-' && gg[g-1] == '-' && ! (gg[g-2] in ['-','.','*'])
                newgg *= '*'
            end
        end
        xx[i,:truncated] = newgg # replace(gg,"-"=> "")

        if occursin(Regex(xx[i,:truncated]),xx[i,3])
            occurstrunc[i] = true
        elseif occursin("N",xx[i,3])
            for letter in ["A","T","C","G"]
                tmp = replace(xx[i,3],"N"=>letter)
                if occursin(Regex(xx[i,:truncated]),tmp)
                    occurstrunc[i] = true
                end
            end
        end

    end
    try
        xx[1,1]
    catch
        # println("fail 4")
        return nothing
    end
    dd = DataFrame(
    PindelID = xx[:,1],
    OccursTrunc = occurstrunc, OccursLiteral = occurslit,
    PindelREF = xx[:,2], PindelALT = xx[:,3],
    ComplexTrunc = xx.truncated,
    ComplexLit = xx.literal, ComplexOld = xx[:,18],
    CallerComplex = xx[:,4], CallerComplexREF= xx[:,5], CallerComplexALT= xx[:,6],
    CallerID = xx[:,7], PindelConfInt = xx[:,8], PindelVAF  = xx[:,9],
    CallerVAFmed = Statistics.median.(xx[:,10]), CallerVAF = sort.(xx[:,10]),
    CHROM  = xx[:,11], SampleTag  = xx[:,12], PindelPOS  = xx[:,13], CallerPOS  = xx[:,14],
    CallerREF= xx[:,15],CallerALT= xx[:,16] ,PindelReadSupport = xx[:,17])
    dd.ID = 1:length(dd.PindelID)
    
    try
        dd[1,1]
    catch
        return nothing
    end
    
    select!(dd, Not([:PindelConfInt,:PindelVAF,:CallerVAF,:CallerVAFmed,:CallerID,:PindelID,:OccursTrunc,:OccursLiteral] ) )
    dd.keeps = trues(nrow(dd))
    for i in 1:nrow(dd)
        l=length(dd[i,"CallerREF"])
        if length(dd[i,"PindelREF"][2:end]) < 2 || length(dd[i,"PindelALT"][2:end]) < 2 || l < 2 
            dd.keeps[i] = false
        elseif dd[i,"CallerPOS"][1] - dd[i,"PindelPOS"] + length(dd[i,"CallerREF"][1]) < 0
            dd.keeps[i]=false
        # elseif 
        #     dd.keeps[i]=false
        else
            pkey = string(dd[i,"PindelPOS"]," ",dd[i,"PindelREF"][2:end]," ",dd[i,"PindelALT"][2:end])
            cpos = dd[i,"CallerPOS"][1]
            for j in 1:l 
                if dd[i,"CallerPOS"][j] < cpos
                    dd.keeps[i]=false
                    break
                end
                cpos += length(dd[i,"CallerREF"][j])
                ckey = string(dd[i,"CallerPOS"][j]," ",dd[i,"CallerREF"][j]," ",dd[i,"CallerALT"][j])
                ckey2 = string(dd[i,"CallerPOS"][j]-1," ",dd[i,"CallerREF"][j]," ",dd[i,"CallerALT"][j])
                ckey3 = string(dd[i,"CallerPOS"][j]+1," ",dd[i,"CallerREF"][j]," ",dd[i,"CallerALT"][j])
                if pkey == ckey || (pkey == ckey2 && (length(dd[i,"CallerREF"][j])>3 && length(dd[i,"CallerALT"][j])>3)) || (pkey == ckey3 && (length(dd[i,"CallerREF"][j])>3 && length(dd[i,"CallerALT"][j])>3))
                    dd.keeps[i]=false
                    break
                end
            end
        end
    end
    filter!(x->x.keeps,dd)
    # dd = dd[sortperm(dd.PindelReadSupport,rev=true),:]
    ddd=unique(select(dd, Not([:ID,:keeps] ) ))
    nrow(ddd) == 0 && return nothing

    return ddd
end

function submain(df,p,oname,P,directoryloc)
    xx=submain2(groupsfun(copy(df),p),P)
    if isnothing(xx)
        CSV.write(directoryloc*"output_"*oname*".tsv.gz",DataFrame(key=["missing","missing"],
        VAF=["missing","missing"],calpos_best=["missing","missing"],calref_best=["missing","missing"],calalt_best=["missing","missing"],bestVAF=["missing","missing"]
        ),delim="\t",compress=true)
        # CSV.write(ARGS[6]*"output_"*oname*".tsv.gz",DataFrame(x1=[0,0],x2=[0,0]),delim="\t",compress=true)
        return nothing
    end
    filter!(x->length(x.calpos_best)==length(x.calref_best)==length(x.calalt_best),xx)
    # println("almost")
    temp2 = df[:,[:key,:emp_vaf]]
    vaf = Dict()
    for i in 1:nrow(temp2)
        # tt = split(temp2[i,:gt_AD],",")
        vaf[temp2[i,:key]] = temp2[i,:emp_vaf] # parse(Float64,string(tt[2])) / (parse(Float64,string(tt[2]))+parse(Float64,string(tt[1])))
    end
    # println(vaf["10979_6_1 chr2 25246632 C>T"])
    xx.bestsamplekeys = [[] for i in 1:nrow(xx)]
    xx.bestVAF = [[] for i in 1:nrow(xx)]
    for i in 1:nrow(xx)
        s=split(xx.CHROM[i],"chr")
        temp = []
        for j in 1:length(xx.calpos_best[i])
            push!(temp,string("chr",s[2]," ",xx.calpos_best[i][j]," ",xx.calref_best[i][j],">",xx.calalt_best[i][j]) )
        end
        xx.bestsamplekeys[i] = temp
        temp = []
        for j in 1:length(xx.bestsamplekeys[i])
            push!(temp, vaf[xx.bestsamplekeys[i][j]] )
        end
        xx.bestVAF[i] = temp
    end
    select!(xx,Not(:bestsamplekeys))
    # println(xx.bestVAF)
    x2 = select(xx,[:key,:VAF,:calpos_best,:calref_best,:calalt_best,:bestVAF])
    CSV.write(directoryloc*"output_"*oname*".tsv.gz",x2,delim="\t",compress=true)
    # CSV.write(ARGS[6]*"output_"*oname*".tsv.gz",x2,delim="\t",compress=true)
    return x2
end

function submain2(d,P)
    temp = []
    if isnothing(d) || nrow(d) == 0
        return nothing
    end

    for g in groupby(d, [:PindelREF, :PindelALT, :CHROM, :SampleTag, :PindelPOS])
        temppos = []
        tempref = []
        tempalt = []
        for i in 1:nrow(g)
            push!(temppos, g[i,:CallerPOS])
            push!(tempalt, g[i,:CallerALT])
            tempalt = replace.(tempalt,"\""=>"")
            push!(tempref, g[i,:CallerREF])
            tempref = replace.(tempref,"\""=>"")
        end
            g = DataFrame(g)
            g.temppos = temppos
            g.tempref = tempref
            g.tempalt = tempalt
            sort!(g,[order(:CallerPOS, by=length, rev=true), order(:PindelReadSupport, rev=true)  ])
            # println(g)
            # println(g[1,:])
            # println(DataFrame(g[1,:]))
            # println(asdfasdf)
        # if length(temppos) < 30
            push!(temp, 
            # DataFrame(pinref=g.PindelREF[1], pinalt=g.PindelALT[1], chrom=g.CHROM[1], sample=g.SampleTag[1], pos=g.PindelPOS[1], calpos=[temppos], calref=[tempref], calalt=[tempalt]) 
            DataFrame(g[1,["PindelREF","PindelALT","CHROM","PindelPOS","temppos","tempref","tempalt"]]))
        # end
    end
    temp = vcat(temp...)
    temp.PindelPOS .+= 1
    temp.PindelREF = chop.(temp.PindelREF,head=1,tail=0)
    temp.PindelALT = chop.(temp.PindelALT,head=1,tail=0)
    # P = CSV.read(qname,DataFrame,header=0; stringtype=String)
    # rename!(P,:Column2 => :CHROM)
    # rename!(P,:Column3 => :POS)
    # rename!(P,:Column4 => :REF)
    # rename!(P,:Column5 => :ALT)
    # rename!(P,:Column6 => :FILTER)
    # rename!(P,:Column7 => :AD)
    # altref = split.(P.gt_AD,",")
    altcount = copy(P.gt_AD_alt_vardict)
    refcount = copy(P.gt_AD_ref_vardict)
    P.VAF = convert.(Float64,altcount) ./ ( convert.(Float64,altcount) .+ convert.(Float64,refcount) )
    P = leftjoin(P, temp, on=[:CHROM => :CHROM, :POS => :PindelPOS, :REF => :PindelREF, :ALT => :PindelALT] )
    # P.altcount = [parse(Int,i[2]) for i in altref]
    # P.refcount = [parse(Int,i[1]) for i in altref]
    # P.binary = [ismissing(i) ? 0 : 1 for i in P.temppos]
    rename!(P, :temppos=>:calpos_best)
    rename!(P, :tempref=>:calref_best)
    rename!(P, :tempalt=>:calalt_best)
    # P.calpos_best = [ismissing(P[i,:temppos]) ? missing : P[i,:temppos][[ ismissing(i) ? missing : argmax(i) for i in [ ismissing(i) ? missing : [sum(length.(j)) for j in i] for i in P[:,:temppos] ]][i]] for i in 1:nrow(P)]
    # P.calref_best = [ismissing(P[i,:tempref]) ? missing : P[i,:tempref][[ ismissing(i) ? missing : argmax(i) for i in [ ismissing(i) ? missing : [sum(length.(j)) for j in i] for i in P[:,:tempref] ]][i]] for i in 1:nrow(P)]
    # P.calalt_best = [ismissing(P[i,:tempalt]) ? missing : P[i,:tempalt][[ ismissing(i) ? missing : argmax(i) for i in [ ismissing(i) ? missing : [sum(length.(j)) for j in i] for i in P[:,:tempalt] ]][i]] for i in 1:nrow(P)]
    return filter(x->!ismissing(x.calref_best),P)
    # return filter(x->length(x.calref_b,P)
    # return P
end

function sqo_main(M,P,oname,directoryloc)
    return submain(subject(copy(M)), query(copy(P)),oname, copy(P),directoryloc)
end


################################################################ CLEAN #################################
function clean(M,L,V,SAMPLE_NAME)
    filter!(x->!occursin("IRC",x.FILTER_Vardict),V);
    filter!(x->!occursin("IRC",x.FILTER_Mutect),M);
    filter!(x->!occursin("IRC",x.FILTER_Lofreq),L);

    


    ## LOFREQ 

    # println(L[1,1:6])
    select!(L,["sample_key","REF","ALT",
    # "VariantClass",
    "QUAL_Lofreq","FILTER_Lofreq",
            "gt_AD_alt_Lofreq","gt_AD_ref_Lofreq","gt_AF_Lofreq","RefFwd","RefRev","AltFwd","AltRev",
            "dust_score","dust_score_3","dust_score_5","dust_score_10","pass_homopolymer_filter",
            # "total_lofreq",
            "FPpass",
            # "total_greater_than_min_alt_count_lofreq",
            "FP_Filter",
            # "total_greater_than_min_VAF_Lofreq",
            "max_gnomAD_AF",
            # "PON_2AT2_percent",
            # "PON_NAT2_percent","PON_MAX_VAF",
            "PON_RefDepth","PON_AltDepth","pon_pvalue","clinvar_SCORE_VEP",
            # "n_samples_Lofreq","n_samples_greater_than_min_VAF_Lofreq","n_samples_greater_than_min_alt_count_lofreq",
            "source.totals.loci","source.totals.p","source.totals.c","CosmicCount","heme_cosmic_count","myeloid_cosmic_count","IsComplex"
            ]);
    rename!(L,:REF=>:REF_combine);
    rename!(L,:ALT=>:ALT_combine);
    rename!(L,:dust_score=>:dust_score_combine);
    rename!(L,:dust_score_3=>:dust_score_3_combine);
    rename!(L,:dust_score_5=>:dust_score_5_combine);
    rename!(L,:dust_score_10=>:dust_score_10_combine);
    rename!(L,:pass_homopolymer_filter=>:pass_homopolymer_filter_combine);
    rename!(L,:max_gnomAD_AF=>:max_gnomAD_AF_combine);
    # rename!(L,:PON_2AT2_percent=>:PON_2AT2_percent_combine);
    # rename!(L,:PON_NAT2_percent=>:PON_NAT2_percent_combine);
    # rename!(L,:PON_MAX_VAF=>:PON_MAX_VAF_combine);
    rename!(L,:PON_RefDepth=>:PON_RefDepth_combine);
    rename!(L,:PON_AltDepth=>:PON_AltDepth_combine);
    rename!(L,:clinvar_SCORE_VEP=>:clinvar_SCORE_VEP_combine);
    # rename!(L,:total_Lofreq=>:n_samples_combine);
    # rename!(L,:total_greater_than_min_VAF_Lofreq=>:n_samples_greater_than_min_VAF_combine);
    # rename!(L,:total_greater_than_min_alt_count_Lofreq=>:n_samples_greater_than_min_alt_count_combine);
    rename!(L,"source.totals.loci"=>"sourcetotalsloci_combine");
    rename!(L,"source.totals.p"=>"sourcetotalsp_combine");
    rename!(L,"source.totals.c"=>"sourcetotalsc_combine");
    rename!(L,"CosmicCount"=>"CosmicCount_combine");
    rename!(L,"heme_cosmic_count"=>"hemecosmiccount_combine");
    rename!(L,"myeloid_cosmic_count"=>"myeloidcosmiccount_combine");
    rename!(L,"IsComplex"=>"IsComplex_combine");



    for nam in names(L)
        if ! occursin("ofreq",nam) && nam != "sample_key" 
            rename!(L,"$(nam)"=>string(nam,"_Lofreq"))
        end
    end

    
    ## Vardict

    select!(V,["sample_key","REF","ALT",
    # "VariantClass",
    "QUAL_Vardict","FILTER_Vardict",
            "PMEAN",
    #         "ReadQual",
            "SBF","ODDRATIO","SN","HIAF","ADJAF","SHIFT3","NM","HICNT","HICOV","gt_AD_alt_Vardict","gt_AD_ref_Vardict",
            "gt_AF_Vardict","RefFwd","RefRev","AltFwd","AltRev","pon_pvalue",
            "dust_score","dust_score_3","dust_score_5","dust_score_10","pass_homopolymer_filter",
            # "total_Vardict",
            "FPpass",
            # "total_greater_than_min_alt_count_Vardict",
            "FP_Filter",
            # "total_greater_than_min_VAF_Vardict",
            "max_gnomAD_AF",
            # "PON_2AT2_percent",
            # "PON_NAT2_percent","PON_MAX_VAF",
            "PON_RefDepth","PON_AltDepth","clinvar_SCORE_VEP",
            # "n_samples_Vardict",
            # "n_samples_greater_than_min_VAF_Vardict","n_samples_greater_than_min_alt_count_Vardict",
            "source.totals.loci","source.totals.p","source.totals.c","CosmicCount","heme_cosmic_count","myeloid_cosmic_count","IsComplex"
            ]);

    rename!(V,:REF=>:REF_combine);
    rename!(V,:ALT=>:ALT_combine);
    # rename!(V,:VariantClass=>:VARIANT_CLASS_combine);
    rename!(V,:dust_score=>:dust_score_combine);
    rename!(V,:dust_score_3=>:dust_score_3_combine);
    rename!(V,:dust_score_5=>:dust_score_5_combine);
    rename!(V,:dust_score_10=>:dust_score_10_combine);
    rename!(V,:pass_homopolymer_filter=>:pass_homopolymer_filter_combine);
    rename!(V,:max_gnomAD_AF=>:max_gnomAD_AF_combine);
    # rename!(V,:PON_2AT2_percent=>:PON_2AT2_percent_combine);
    # rename!(V,:PON_NAT2_percent=>:PON_NAT2_percent_combine);
    # rename!(V,:PON_MAX_VAF=>:PON_MAX_VAF_combine);
    rename!(V,:PON_RefDepth=>:PON_RefDepth_combine);
    rename!(V,:PON_AltDepth=>:PON_AltDepth_combine);
    rename!(V,:clinvar_SCORE_VEP=>:clinvar_SCORE_VEP_combine);
    # rename!(V,:total_Vardict=>:n_samples_combine);
    # rename!(V,:total_greater_than_min_VAF_Vardict=>:n_samples_greater_than_min_VAF_combine);
    # rename!(V,:total_greater_than_min_alt_count_Vardict=>:n_samples_greater_than_min_alt_count_combine);
    rename!(V,"source.totals.loci"=>"sourcetotalsloci_combine");
    rename!(V,"source.totals.p"=>"sourcetotalsp_combine");
    rename!(V,"source.totals.c"=>"sourcetotalsc_combine");
    rename!(V,"CosmicCount"=>"CosmicCount_combine");
    rename!(V,"heme_cosmic_count"=>"hemecosmiccount_combine");
    rename!(V,"myeloid_cosmic_count"=>"myeloidcosmiccount_combine");
    rename!(V,"IsComplex"=>"IsComplex_combine");

    for nam in names(V)
        if ! occursin("ardict",nam) && nam != "sample_key" 
            rename!(V,"$(nam)"=>string(nam,"_Vardict")) 
        end
    end

    ## MUTECTOR

    select!(M,["sample_key","REF","ALT",
    # "VariantClass",
    # "QUAL_Mutect",
    "FILTER_Mutect",
            "GERMQ","MBQ","MFRL","MMQ","MPOS","ROQ","TLOD","gt_AD_alt_Mutect","gt_AD_ref_Mutect",
            "gt_AF_Mutect","RefFwd","RefRev","AltFwd","AltRev","pon_pvalue",
            # "total_greater_than_min_alt_count_Mutect",
            "dust_score","dust_score_3","dust_score_5","dust_score_10","pass_homopolymer_filter",
            # "total_Mutect",
            "FPpass","FP_Filter"  ,
            # "total_greater_than_min_VAF_Mutect",
            "max_gnomAD_AF",
            # "PON_2AT2_percent",
            # "PON_NAT2_percent","PON_MAX_VAF",
            "PON_RefDepth","PON_AltDepth","clinvar_SCORE_VEP",
            # "n_samples_Mutect",
            # "n_samples_greater_than_min_VAF_Mutect","n_samples_greater_than_min_alt_count_Mutect",
            "source.totals.loci","source.totals.p","source.totals.c","CosmicCount","heme_cosmic_count","myeloid_cosmic_count","IsComplex"
            ]);

    rename!(M,:REF=>:REF_combine);
    rename!(M,:ALT=>:ALT_combine);
    # rename!(M,:VariantClass=>:VARIANT_CLASS_combine);
    rename!(M,:dust_score=>:dust_score_combine);
    rename!(M,:dust_score_3=>:dust_score_3_combine);
    rename!(M,:dust_score_5=>:dust_score_5_combine);
    rename!(M,:dust_score_10=>:dust_score_10_combine);
    rename!(M,:pass_homopolymer_filter=>:pass_homopolymer_filter_combine);
    rename!(M,:max_gnomAD_AF=>:max_gnomAD_AF_combine);
    # rename!(M,:PON_2AT2_percent=>:PON_2AT2_percent_combine);
    # rename!(M,:PON_NAT2_percent=>:PON_NAT2_percent_combine);
    # rename!(M,:PON_MAX_VAF=>:PON_MAX_VAF_combine);
    rename!(M,:PON_RefDepth=>:PON_RefDepth_combine);
    rename!(M,:PON_AltDepth=>:PON_AltDepth_combine);
    rename!(M,:clinvar_SCORE_VEP=>:clinvar_SCORE_VEP_combine);
    # rename!(M,:total_Mutect=>:n_samples_combine);
    # rename!(M,:total_greater_than_min_VAF_Mutect=>:n_samples_greater_than_min_VAF_combine);
    # rename!(M,:total_greater_than_min_alt_count_Mutect=>:n_samples_greater_than_min_alt_count_combine);
    rename!(M,"source.totals.loci"=>"sourcetotalsloci_combine");
    rename!(M,"source.totals.p"=>"sourcetotalsp_combine");
    rename!(M,"source.totals.c"=>"sourcetotalsc_combine");
    rename!(M,"CosmicCount"=>"CosmicCount_combine");
    rename!(M,"heme_cosmic_count"=>"hemecosmiccount_combine");
    rename!(M,"myeloid_cosmic_count"=>"myeloidcosmiccount_combine");
    rename!(M,"IsComplex"=>"IsComplex_combine");

    for nam in names(M)
        if ! occursin("utect",nam) && nam != "sample_key" 
            rename!(M,"$(nam)"=>string(nam,"_Mutect"))
        end
    end

    ## COMBINE

    temp = outerjoin(V,L,on=[:sample_key])
    temp = outerjoin(temp,M,on=[:sample_key]);

    temp.IsComplex_combine_Mutect = [ismissing(i) ? false : i for i in temp.IsComplex_combine_Mutect]
    temp.IsComplex_combine_Lofreq = [ismissing(i) ? false : i for i in temp.IsComplex_combine_Lofreq]
    temp.IsComplex_combine_Vardict = [ismissing(i) ? false : i for i in temp.IsComplex_combine_Vardict]
    combiner!(temp,"IsComplex_combine")
    combiner!(temp,"REF_combine")
    combiner!(temp,"ALT_combine")
    # temp = combiner!(temp,"VARIANT_CLASS_combine")
    combiner!(temp,"dust_score_combine")
    combiner!(temp,"dust_score_3_combine")
    combiner!(temp,"dust_score_5_combine")
    combiner!(temp,"dust_score_10_combine")
    combiner!(temp,"pass_homopolymer_filter_combine")
    combiner!(temp,"max_gnomAD_AF_combine")
    # temp = combiner!(temp,"PON_2AT2_percent_combine")
    # combiner!(temp,"PON_NAT2_percent_combine")
    # combiner!(temp,"PON_MAX_VAF_combine")
    combiner!(temp,"PON_RefDepth_combine")
    combiner!(temp,"PON_AltDepth_combine")
    combiner!(temp,"clinvar_SCORE_VEP_combine")
    # temp = combiner!(temp,"n_samples_combine")
    # temp = combiner!(temp,"n_samples_greater_than_min_VAF_combine")
    # temp = combiner!(temp,"n_samples_greater_than_min_alt_count_combine")
    combiner!(temp,"sourcetotalsloci_combine")
    combiner!(temp,"sourcetotalsp_combine")
    combiner!(temp,"sourcetotalsc_combine")
    combiner!(temp,"CosmicCount_combine")
    combiner!(temp,"hemecosmiccount_combine")
    combiner!(temp,"myeloidcosmiccount_combine")



    # temp.PON2AT2percent = [ismissing(i) ? 0 : i for i in temp.PON2AT2percent];
    # temp.PONNAT2percent = [ismissing(i) ? 0 : i for i in temp.PONNAT2percent];
    # temp.PONMAXVAF = [ismissing(i) ? missing : i for i in temp.PONMAXVAF];
    temp.PONRefDepth = [ismissing(i) ? missing : i for i in temp.PONRefDepth];
    temp.PONAltDepth = [ismissing(i) ? missing : i for i in temp.PONAltDepth];
    temp.clinvarSCOREVEP = [ismissing(i) ? 0 : i for i in temp.clinvarSCOREVEP];
    # println(temp.clinvarSCOREVEP)
    # temp.nsamples = [ismissing(i) ? 1 : i for i in temp.nsamples];
    # temp.nsamplesgreaterthanminVAF = [ismissing(i) ? 0 : i for i in temp.nsamplesgreaterthanminVAF];
    # temp.nsamplesgreaterthanminaltcount = [ismissing(i) ? 0 : i for i in temp.nsamplesgreaterthanminaltcount];
    temp.sourcetotalsloci = [ismissing(i) ? 0 : i for i in temp.sourcetotalsloci];
    temp.sourcetotalsp = [ismissing(i) ? 0 : i for i in temp.sourcetotalsp];
    temp.sourcetotalsc = [ismissing(i) ? 0 : i for i in temp.sourcetotalsc];
    temp.CosmicCount = [ismissing(i) ? 0 : i for i in temp.CosmicCount];
    temp.hemecosmiccount = [ismissing(i) ? 0 : i for i in temp.hemecosmiccount];
    temp.myeloidcosmiccount = [ismissing(i) ? 0 : i for i in temp.myeloidcosmiccount];


    temp.ref_len = float.(length.(temp.REF))
    temp.alt_len = float.(length.(temp.ALT))

    if nrow(M)>0
        temp.gt_AD_ref_Mutect_z = abs.(temp.gt_AD_ref_Mutect .- Statistics.mean(skipmissing(temp.gt_AD_ref_Mutect))) ./ Statistics.std(skipmissing(temp.gt_AD_ref_Mutect));
    else
        temp.gt_AD_ref_Mutect_z = [missing for i in 1:nrow(temp)]
    end
    temp.gt_AD_ref_Vardict_z = abs.(temp.gt_AD_ref_Vardict .- Statistics.mean(skipmissing(temp.gt_AD_ref_Vardict))) ./ Statistics.std(skipmissing(temp.gt_AD_ref_Vardict));
    if nrow(L)>0
        temp.gt_AD_ref_Lofreq_z = abs.(temp.gt_AD_ref_Lofreq .- Statistics.mean(skipmissing(temp.gt_AD_ref_Lofreq))) ./ Statistics.std(skipmissing(temp.gt_AD_ref_Lofreq));
    else
        temp.FILTER_Lofreq_PASS = [missing for i in 1:nrow(temp)]
        temp.gt_AD_ref_Lofreq_z = [missing for i in 1:nrow(temp)]
    end

    SBF_pvalue_Mutect = zeros(nrow(temp))
    for i in 1:nrow(temp)
        try
            SBF_pvalue_Mutect[i] = HypothesisTests.pvalue(HypothesisTests.FisherExactTest(
                    temp.RefFwd_Mutect[i],
                    temp.RefRev_Mutect[i],
                    temp.AltFwd_Mutect[i],
                    temp.AltRev_Mutect[i]
                    ))
        catch 
            SBF_pvalue_Mutect[i] = NaN
        end
    end
    temp.SBF_pvalue_Mutect = SBF_pvalue_Mutect

    SBF_pvalue_Vardict = zeros(nrow(temp))
    for i in 1:nrow(temp)
        try
            SBF_pvalue_Vardict[i] = HypothesisTests.pvalue(HypothesisTests.FisherExactTest(
                    temp.RefFwd_Vardict[i],
                    temp.RefRev_Vardict[i],
                    temp.AltFwd_Vardict[i],
                    temp.AltRev_Vardict[i]
                    ))
        catch 
            SBF_pvalue_Vardict[i] = NaN
        end
    end
    temp.SBF_pvalue_Vardict = SBF_pvalue_Vardict

    SBF_pvalue_Lofreq = zeros(nrow(temp))
    for i in 1:nrow(temp)
        try
            SBF_pvalue_Lofreq[i] = HypothesisTests.pvalue(HypothesisTests.FisherExactTest(
                    temp.RefFwd_Lofreq[i],
                    temp.RefRev_Lofreq[i],
                    temp.AltFwd_Lofreq[i],
                    temp.AltRev_Lofreq[i]
                    ))
        catch 
            SBF_pvalue_Lofreq[i] = NaN
        end
    end
    temp.SBF_pvalue_Lofreq = SBF_pvalue_Lofreq

    temp.REF = [length(i)>1 ? "I" : i for i in temp.REF];
    temp.ALT = [length(i)>1 ? "I" : i for i in temp.ALT];


    ## V

    # v = copy(temp);
    temp.passhomopolymerfilter = [ismissing(i) ? missing : i for i in temp.passhomopolymerfilter]
    temp.FPpass_Mutect = [ismissing(i) ? missing : i for i in temp.FPpass_Mutect]
    temp.FPpass_Lofreq = [ismissing(i) ? missing : i for i in temp.FPpass_Lofreq]
    temp.FPpass_Vardict = [ismissing(i) ? missing : i for i in temp.FPpass_Vardict]
    # temp.real = [ismissing(i) ? missing : i for i in temp.real]

    temp.QUAL_Vardict = float.(temp.QUAL_Vardict)
    temp.PMEAN_Vardict = float.(temp.PMEAN_Vardict)
    # temp.ReadQual_Vardict = float.(temp.ReadQual_Vardict)
    temp.SBF_Vardict = float.(temp.SBF_Vardict)
    temp.ODDRATIO_Vardict = float.(temp.ODDRATIO_Vardict)
    temp.SN_Vardict = float.(temp.SN_Vardict)
    temp.HIAF_Vardict = float.(temp.HIAF_Vardict)
    temp.ADJAF_Vardict = float.(temp.ADJAF_Vardict)
    temp.SHIFT3_Vardict = float.(temp.SHIFT3_Vardict)
    temp.NM_Vardict = float.(temp.NM_Vardict)
    temp.HICNT_Vardict = float.(temp.HICNT_Vardict)
    temp.HICOV_Vardict = float.(temp.HICOV_Vardict)
    temp.gt_AD_ref_Vardict = float.(temp.gt_AD_ref_Vardict)
    temp.gt_AD_alt_Vardict = float.(temp.gt_AD_alt_Vardict)
    temp.gt_AF_Vardict = float.(temp.gt_AF_Vardict)
    # temp.pon_pvalue_Vardict = float.(temp.pon_pvalue_Vardict)
    temp.RefFwd_Vardict = float.(temp.RefFwd_Vardict)
    temp.RefRev_Vardict = float.(temp.RefRev_Vardict)
    temp.AltFwd_Vardict = float.(temp.AltFwd_Vardict)
    temp.AltRev_Vardict = float.(temp.AltRev_Vardict)


    temp.QUAL_Lofreq = float.(temp.QUAL_Lofreq)
    temp.gt_AD_ref_Lofreq = float.(temp.gt_AD_ref_Lofreq)
    temp.gt_AD_alt_Lofreq = float.(temp.gt_AD_alt_Lofreq)
    temp.gt_AF_Lofreq = float.(temp.gt_AF_Lofreq)
    # temp.pon_pvalue_Lofreq = float.(temp.pon_pvalue_Lofreq)
    temp.RefFwd_Lofreq = float.(temp.RefFwd_Lofreq)
    temp.RefRev_Lofreq = float.(temp.RefRev_Lofreq)
    temp.AltFwd_Lofreq = float.(temp.AltFwd_Lofreq)
    temp.AltRev_Lofreq = float.(temp.AltRev_Lofreq)
            
    temp.GERMQ_Mutect = float.(temp.GERMQ_Mutect)
    temp2 = [ismissing(i) ? [missing,missing] : split(i,",") for i in temp.MBQ_Mutect]
    temp.MBQ_Mutect_1 = [ismissing(i[1]) ? missing : parse(Float64,i[1]) for i in temp2]
    temp.MBQ_Mutect_2 = [ismissing(i[2]) ? missing : parse(Float64,i[2]) for i in temp2]

    temp2 = [ismissing(i) ? [missing,missing] : split(i,",") for i in temp.MFRL_Mutect]
    temp.MFRL_Mutect_1 = [ismissing(i[1]) ? missing : parse(Float64,i[1]) for i in temp2]
    temp.MFRL_Mutect_2 = [ismissing(i[2]) ? missing : parse(Float64,i[2]) for i in temp2]

    temp2 = [ismissing(i) ? [missing,missing] : split(i,",") for i in temp.MMQ_Mutect]
    temp.MMQ_Mutect_1 = [ismissing(i[1]) ? missing : parse(Float64,i[1]) for i in temp2]
    temp.MMQ_Mutect_2 = [ismissing(i[2]) ? missing : parse(Float64,i[2]) for i in temp2]

    temp.MPOS_Mutect = float.(temp.MPOS_Mutect)
    temp.ROQ_Mutect = float.(temp.ROQ_Mutect)
    temp.TLOD_Mutect = float.(temp.TLOD_Mutect)
    temp.gt_AD_ref_Mutect = float.(temp.gt_AD_ref_Mutect)
    temp.gt_AD_alt_Mutect = float.(temp.gt_AD_alt_Mutect)
    temp.gt_AF_Mutect = float.(temp.gt_AF_Mutect)
    # temp.pon_pvalue_Mutect = float.(temp.pon_pvalue_Mutect)
    temp.RefFwd_Mutect = float.(temp.RefFwd_Mutect)
    temp.RefRev_Mutect = float.(temp.RefRev_Mutect)
    temp.AltFwd_Mutect = float.(temp.AltFwd_Mutect);
    temp.AltRev_Mutect = float.(temp.AltRev_Mutect);

    temp.gt_AD_ref_Mutect_z = float.(temp.gt_AD_ref_Mutect_z)
    temp.gt_AD_ref_Vardict_z  = float.(temp.gt_AD_ref_Vardict_z)
    temp.gt_AD_ref_Lofreq_z  = float.(temp.gt_AD_ref_Lofreq_z)

    # temp.gt_AD_alt_Mutect_z  = float.(temp.gt_AD_alt_Mutect_z)
    # temp.gt_AD_alt_Vardict_z  = float.(temp.gt_AD_alt_Vardict_z)
    # temp.gt_AD_alt_Lofreq_z = float.(temp.gt_AD_alt_Lofreq_z)

    select!(temp, Not([:MBQ_Mutect,:MFRL_Mutect,:MMQ_Mutect]));




    temp.FILTER_Lofreq = [ismissing(i) ? "Not_Detected" : i=="Not Detected" ? "Not_Detected" : i for i in temp.FILTER_Lofreq];
    temp.FILTER_Mutect = [ismissing(i) ? "Not_Detected" : i=="Not Detected" ? "Not_Detected" : i for i in temp.FILTER_Mutect];
    temp.FILTER_Vardict = [ismissing(i) ? "Not_Detected" : i=="Not Detected" ? "Not_Detected" : i for i in temp.FILTER_Vardict];

    temp.FP_Filter_Lofreq = [ismissing(i) ? "Not_Detected" : i=="Not Detected" ? "Not_Detected" : i for i in temp.FP_Filter_Lofreq];
    temp.FP_Filter_Mutect = [ismissing(i) ? "Not_Detected" : i=="Not Detected" ? "Not_Detected" : i for i in temp.FP_Filter_Mutect];
    temp.FP_Filter_Vardict = [ismissing(i) ? "Not_Detected" : i=="Not Detected" ? "Not_Detected" : i for i in temp.FP_Filter_Vardict];

    temp.FPpass = ones(Float32,nrow(temp));
    temp.FP_Filter = fill("PASS",nrow(temp));


    for i in 1:nrow(temp)
        if temp.FP_Filter_Mutect[i] != "PASS" && temp.FP_Filter_Mutect[i] != "Not_Detected"
            temp.FP_Filter[i] = temp.FP_Filter_Mutect[i]
            temp.FPpass[i] = 0.0f0
        elseif temp.FP_Filter_Lofreq[i] != "PASS" && temp.FP_Filter_Lofreq[i] != "Not_Detected"
            temp.FP_Filter[i] = temp.FP_Filter_Lofreq[i]
            temp.FPpass[i] = 0.0f0
        elseif temp.FP_Filter_Vardict[i] != "PASS" && temp.FP_Filter_Vardict[i] != "Not_Detected"
            temp.FP_Filter[i] = temp.FP_Filter_Vardict[i]
            temp.FPpass[i] = 0.0f0
        end
    end

    select!(temp, Not([
                "FPpass_Mutect",
                "FPpass_Lofreq",
                "FPpass_Vardict",
                "FP_Filter_Mutect",
                "FP_Filter_Vardict",
                "FP_Filter_Lofreq",            
                ]));

    temp.AF = copy(temp.gt_AF_Vardict);
    for i in 1:nrow(temp)
        if ismissing(temp.AF[i])
            if ! ismissing(temp[i,"gt_AF_Lofreq"])
                temp.AF[i] = temp[i,"gt_AF_Lofreq"]
            elseif ! ismissing(temp[i,"gt_AF_Mutect"])
                temp.AF[i] = temp[i,"gt_AF_Mutect"]
            end
        end
    end   

    temp.altdepth = copy(temp.gt_AD_alt_Vardict);
    for i in 1:nrow(temp)
        if ismissing(temp.altdepth[i])
            if ! ismissing(temp[i,"gt_AD_alt_Lofreq"])
                temp.altdepth[i] = temp[i,"gt_AD_alt_Lofreq"]
            elseif ! ismissing(temp[i,"gt_AD_alt_Mutect"])
                temp.altdepth[i] = temp[i,"gt_AD_alt_Mutect"]
            end
        end
    end   

    # v3 = copy(v);


    temp.key = [join(i[2:end]," ") for i in split.(temp.sample_key," ")]
    temp = leftjoin(temp,getpileup(SAMPLE_NAME),on=[:key])

    temp.pon_af_zscore = zeros(Float64,nrow(temp))
    for i in 1:nrow(temp)
        # if ismissing(temp.PONvafstd[i])
        #     println(temp.sample_key[i])
        #     println(temp.key[i])
        #     println(temp.FILTER_Lofreq[i])
        #     println(temp.FILTER_Vardict[i])
        #     println(temp.FILTER_Mutect[i])
        #     println(temp.ALT[i])
        #     println(temp.REF[i])
        #     println(temp.FP_Filter[i])

        # end
        std = temp.PONvafstd[i] < 0.000001 ? 0.000001 : temp.PONvafstd[i]
       temp.pon_af_zscore[i] = (temp.AF[i] - temp.PONvafmean[i]) / std
    end

    select!(temp,Not([
                
    "AltFwd_Mutect",
    "AltRev_Mutect",
    "AltFwd_Lofreq",
    "AltRev_Lofreq",
    "AltFwd_Vardict" ,
    "AltRev_Vardict",
                
    "RefFwd_Mutect",
    "RefRev_Mutect",
    "RefFwd_Lofreq",
    "RefRev_Lofreq",
    "RefFwd_Vardict" ,
    "RefRev_Vardict",
                
                "SBF_Vardict",
                "gt_AD_alt_Mutect",
                "gt_AD_alt_Vardict",
                "gt_AD_alt_Lofreq",
                
                "gt_AF_Mutect",
                "gt_AF_Vardict",
                "gt_AF_Lofreq",
                
                "ADJAF_Vardict",
                "HIAF_Vardict",
                # "AF",
                # "altdepth","key","sample",
                # "pon_mutect","pon_lofreq","pon_vardict","totaldepth_sum","altdepth_sum","vaf_std","vaf_mean",
                # "PON_2AT2_percent",
                ]) );


    # coldrop = []

    # println("dropping")
    # for i in names(v4)
    # #     if i == "col1" || i=="training"
    # #         continue
    # #     end
    #     if length(unique(v4[:,i])) < 2 || occursin("Pindel",i)
    #         push!(coldrop,i)
    #         println(i)
    #     end 
    # end
    # unique!(coldrop)
    # select!(v4, Not(coldrop));
    # x=copy(v4);



    for j in 1:nrow(temp)
        if typeof(temp.sourcetotalsloci[j]) == String 
            temp.sourcetotalsloci[j] = sum([parse(Float32,i[2]) for i in split.( split(temp.sourcetotalsloci[j],","), ":")])
        end
        if typeof(temp.sourcetotalsp[j]) == String 
            temp.sourcetotalsp[j] = sum([parse(Float32,i[2]) for i in split.( split(temp.sourcetotalsp[j],","), ":")])
        end
        if typeof(temp.sourcetotalsc[j]) == String 
            temp.sourcetotalsc[j] = sum([parse(Float32,i[2]) for i in split.( split(temp.sourcetotalsc[j],","), ":")])
        end
    end


    temp = hcat(
        temp,
        textonehot(temp,"FILTER_Mutect"),
        textonehot(temp,"FILTER_Lofreq"),
        textonehot(temp,"FILTER_Vardict"),
        textonehot(temp,"FP_Filter"),
    )

    select!(temp,Not([
                "FILTER_Lofreq","FILTER_Vardict","FILTER_Mutect",
                "ALT","REF",
                # "VARIANTCLASS",
                "FP_Filter",

                # "AF",
                # "altdepth"
                ]) );

    allowmissing!(temp);
    for col in 1:ncol(temp)
        temp[!,col] = [ismissing(i) ? missing : typeof(i) == String ? i : Float32(i) for i in temp[:,col] ]
    end
    for j in 1:ncol(temp)
        for i in 1:nrow(temp)
            if typeof(temp[i,j])!=String && (ismissing(temp[i,j]) || isnan(temp[i,j]) || isinf(temp[i,j]))
                temp[i,j] = NaN32
            end
        end
        if occursin("pvalue",names(temp)[j])
            for i in 1:nrow(temp)
                if ! ismissing(temp[i,j])
                    if isapprox(temp[i,j], 0.0) 
                        temp[i,j] = -307.65
                    else
                        temp[i,j] = log10(temp[i,j])
                    end
                end
            end
        end
    end

    for col in 1:ncol(temp)
        temp[!,col] = [ismissing(i) ? missing : typeof(i) == String ? i : Float32(i) for i in temp[:,col] ]
    end

    return temp

end
################################################################ PON #################################
function getpileup(subject)
    df=CSV.read(subject,DataFrame,comment="##",delim="\t",buffer_in_memory=true,ntasks=1,truestrings=["true", "True", "TRUE", "1"],falsestrings=["false", "False", "FALSE", "0"])
    df.key = ["" for i in 1:nrow(df)]
    for i in 1:nrow(df)
        df.key[i] = string(df[i,1]," ",df[i,2]," ",df[i,4],">",df[i,5])
    end
    df.PONRefDepth = zeros(nrow(df))
    df.PONAltDepth = zeros(nrow(df))
    df.PONvafmean = zeros(nrow(df))
    df.PONvafstd = zeros(nrow(df))
    temp = [0.0 for i in 1:27]
    for i in 1:nrow(df)
        for k in 1:27
            temp[k] = 0.0
        end
        for j in 10:ncol(df)-5
            pp = split(df[i,j],":")
            df.PONRefDepth[i] += parse(Float32,pp[1]) # DP total depths
            df.PONAltDepth[i] += parse(Float32,pp[3]) # ad alt depths
            df.PONvafmean[i] += parse(Float64,pp[4]) # vaf 
            temp[j-9] = parse(Float64,pp[4])
        end
        df.PONvafmean[i] /= 27
        df.PONvafstd[i] = Statistics.std(temp)
    end
    return unique(select(df,[:key,:PONvafmean,:PONvafstd])) # df[:,end-4:end]
end
################################################################ PREDICT #################################
function predict(df,modelloc)
    Random.seed!(Dates.now().instant.periods.value);
    # mach = MLJ.machine(ARGS[7]*"$(rand(0:4))XGB_FINAL_MODEL.jlso")
    mach = MLJ.machine(modelloc)
    sample_key = copy(df.sample_key)
    key = copy(df.key)
    AF = copy(df.AF)
    # altdepth = copy(df.altdepth)
    # x = select(df, Not([:key,:sample_key,:altdepth,:AF]) );
    # finalnames = CSV.read("features_final.csv.gz",DataFrame,ntasks=1,truestrings=["true", "True", "TRUE", "1"],falsestrings=["false", "False", "FALSE", "0"]).name
    # finalnames = CSV.read(ARGS[8],DataFrame,ntasks=1,truestrings=["true", "True", "TRUE", "1"],falsestrings=["false", "False", "FALSE", "0"]).name
    finalnames = ["sample_key",
        "QUAL_Vardict",
        "PMEAN_Vardict",
        "ODDRATIO_Vardict",
        "SN_Vardict",
        "SHIFT3_Vardict",
        "NM_Vardict",
        "HICNT_Vardict",
        "HICOV_Vardict",
        "gt_AD_ref_Vardict",
        "pon_pvalue_Vardict",
        "QUAL_Lofreq",
        "gt_AD_ref_Lofreq",
        "pon_pvalue_Lofreq",
        "GERMQ_Mutect",
        "MPOS_Mutect",
        "ROQ_Mutect",
        "TLOD_Mutect",
        "gt_AD_ref_Mutect",
        "pon_pvalue_Mutect",
        "dustscore",
        "dustscore3",
        "dustscore5",
        "dustscore10",
        "passhomopolymerfilter",
        "maxgnomADAF",
        "PONRefDepth",
        "PONAltDepth",
        "clinvarSCOREVEP",
        "sourcetotalsloci",
        "sourcetotalsp",
        "sourcetotalsc",
        "CosmicCount",
        "hemecosmiccount",
        "myeloidcosmiccount",
        "ref_len",
        "alt_len",
        "gt_AD_ref_Mutect_z",
        "gt_AD_ref_Vardict_z",
        "gt_AD_ref_Lofreq_z",
        "SBF_pvalue_Mutect",
        "SBF_pvalue_Vardict",
        "SBF_pvalue_Lofreq",
        "MBQ_Mutect_1",
        "MBQ_Mutect_2",
        "MFRL_Mutect_1",
        "MFRL_Mutect_2",
        "MMQ_Mutect_1",
        "MMQ_Mutect_2",
        "FPpass",
        "key",
        "PONvafmean",
        "PONvafstd",
        "pon_af_zscore",
        "FILTER_Mutect_Not_Detected",
        "FILTER_Mutect_clustered_events",
        "FILTER_Mutect_haplotype",
        "FILTER_Mutect_orientation",
        "FILTER_Mutect_position",
        "FILTER_Mutect_strand_bias",
        "FILTER_Mutect_PASS",
        "FILTER_Mutect_weak_evidence",
        "FILTER_Lofreq_Not_Detected",
        "FILTER_Lofreq_PASS",
        "FILTER_Vardict_PASS",
        "FILTER_Vardict_NM5_25",
        "FILTER_Vardict_p8",
        "FILTER_Vardict_pSTD",
        "FILTER_Vardict_BCBIO",
        "FILTER_Vardict_Not_Detected",
        "FP_Filter_PASS",
        "FP_Filter_RLD25",
        "FP_Filter_DETP20",
        "FP_Filter_MVC4",
        "FP_Filter_SB1",
        "FP_Filter_MMQSD50",
        "FP_Filter_NRC",
        "FP_Filter_PB10",
        "FP_Filter_MMQS100"
    ]
    missing_names = setdiff(finalnames,names(df)) # columns not in finalnames
    if ! isempty(missing_names)
        # if !("FILTER_Vardict_PASS" in names(df))
        #     df[:,"FILTER_Vardict_PASS"] = [ismissing(i) ? NaN32 : 0.0f0 for i in df[:,"pon_pvalue_Vardict"] ]
        # end
        # if !("FILTER_Mutect_PASS" in names(df))
        #     df[:,"FILTER_Mutect_PASS"] = [ismissing(i) ? NaN32 : 0.0f0 for i in df[:,"pon_pvalue_Mutect"] ]
        # end
        # if !("FILTER_Lofreq_PASS" in names(df))
        #     df[:,"FILTER_Lofreq_PASS"] = [ismissing(i) ? NaN32 : 0.0f0 for i in df[:,"pon_pvalue_Lofreq"] ]
        # end
        # println(missing_names)
        for nam in missing_names
            if occursin("FILTER",nam)
                if occursin("utect",nam)
                    df[:,nam] = [ismissing(i) ? NaN32 : 0.0f0 for i in df[:,"pon_pvalue_Mutect"]]
                elseif occursin("ardict",nam)
                    df[:,nam] = [ismissing(i) ? NaN32 : 0.0f0 for i in df[:,"pon_pvalue_Vardict"]]
                else occursin("ofreq",nam)
                    df[:,nam] = [ismissing(i) ? NaN32 : 0.0f0 for i in df[:,"pon_pvalue_Lofreq"]]
                end
            elseif occursin("FP_F",nam)
                df[:,nam] = [0.0f0 for i in 1:nrow(df)]
            else
                df[:,nam] = [NaN32 for i in 1:nrow(df)]
            end
        end
    end
    # allowmissing!(df);
    x = select(df, filter(x->!occursin("key",x),finalnames) );
    for j in 1:ncol(x)
        for i in 1:nrow(x)
            if ismissing(x[i,j])
                try
                x[i,j] = NaN32
                catch 
                    println(typeof(x[:,j]))
                    println(names(x)[j])
                    println(asdfasd)
                end
            end
        end
    end
    disallowmissing!(x);

    x.predicted_probability = vec(MLJ.pdf(MLJ.predict(mach, x ),[true]));
    
    x.AF = AF
    x.sample_key = sample_key
    x.key = key

    
   for i in names(x)
        if ! occursin("key",i)
            rename!(x,i=>string(i,"_XGB"))
        end
    end

   return x #[:,[:sample_key,:key,:predicted_probability,:FP_probability]]
end
################################################################ MAIN #################################
import CSV
using DataFrames
import Statistics
import Combinatorics
import HypothesisTests
import XGBoost
import MLJ
import Dates
import Random
import StringDistances
import StatsBase
# import Arrow
XGB = MLJ.@load XGBoostClassifier verbosity=0 

function main(lofreqname::String,mutectname::String,vardictname::String,pindelname::String,ponname::String,directoryloc::String,modelloc::String)::Nothing

    SAMPLE_NAME = split(mutectname,".")[2] # "218281_1_1000_S212"

    L=CSV.read(directoryloc*lofreqname,DataFrame,delim="\t",downcast=false,stringtype=String, missingstring=["","NA", "NAN", "NULL"],ntasks=1, escapechar='\\',truestrings=["true", "True", "TRUE", "1"],falsestrings=["false", "False", "FALSE", "0"]);
    M=CSV.read(directoryloc*mutectname,DataFrame,delim="\t",downcast=false,stringtype=String, missingstring=["","NA", "NAN", "NULL"],ntasks=1, escapechar='\\',truestrings=["true", "True", "TRUE", "1"],falsestrings=["false", "False", "FALSE", "0"]);
    V=CSV.read(directoryloc*vardictname,DataFrame,delim="\t",downcast=false,stringtype=String, missingstring=["","NA", "NAN", "NULL"],ntasks=1, escapechar='\\',truestrings=["true", "True", "TRUE", "1"],falsestrings=["false", "False", "FALSE", "0"]);
    P = CSV.read(directoryloc*pindelname,DataFrame,delim="\t",comment="##",header=1,ntasks=1,truestrings=["true", "True", "TRUE", "1"],falsestrings=["false", "False", "FALSE", "0"])
    rename!(P,"$SAMPLE_NAME"=>"VALUES")
    rename!(P,"#CHROM"=>"CHROM");
    P.VALUES = [split(split(i,":")[2],",") for i in P.VALUES];
    P.gt_AD_ref = [parse(Int,i[1]) for i in P.VALUES];
    P.gt_AD_alt = [parse(Int,i[2]) for i in P.VALUES];
    select!(P,[:CHROM,:POS,:REF,:ALT,:gt_AD_ref,:gt_AD_alt]);

    rename!(L,:SAMPLE => :SN_TAG);
    rename!(M,:SAMPLE => :SN_TAG);
    rename!(V,:SAMPLE => :SN_TAG);

    sntag!(L.SN_TAG);
    sntag!(M.SN_TAG);
    sntag!(V.SN_TAG);

    # L.sample_key = [string(L.SN_TAG[i]," ",L.CHROM[i]," " , L.POS[i] , " " , L.REF[i] , ">" , L.ALT[i]) for i in 1:nrow(L)]
    # M.sample_key = [string(M.SN_TAG[i]," ",M.CHROM[i]," " , M.POS[i] , " " , M.REF[i] , ">" , M.ALT[i]) for i in 1:nrow(M)]
    # V.sample_key = [string(V.SN_TAG[i]," ",V.CHROM[i]," " , V.POS[i] , " " , V.REF[i] , ">" , V.ALT[i]) for i in 1:nrow(V)]

    L.key = [string(L.CHROM[i]," " , L.POS[i] , " " , L.REF[i] , ">" , L.ALT[i]) for i in 1:nrow(L)]
    M.key = [string(M.CHROM[i]," " , M.POS[i] , " " , M.REF[i] , ">" , M.ALT[i]) for i in 1:nrow(M)]
    V.key = [string(V.CHROM[i]," " , V.POS[i] , " " , V.REF[i] , ">" , V.ALT[i]) for i in 1:nrow(V)]
    P.key = [string(P.CHROM[i]," " , P.POS[i] , " " , P.REF[i] , ">" , P.ALT[i]) for i in 1:nrow(P)]
   

    #### COMPLEX ####
    M.Column1 = 1:nrow(M)
    L.Column1 = 1:nrow(L)
    V.Column1 = 1:nrow(V)
    V.IsComplex = falses(nrow(V))
    M.IsComplex = falses(nrow(M))
    L.IsComplex = falses(nrow(L))
    V.IsComplexPart = falses(nrow(V))
    M.IsComplexPart = falses(nrow(M))
    L.IsComplexPart = falses(nrow(L))
    V.ComplexKey = ["" for _ in 1:nrow(V)]
    M.ComplexKey = ["" for _ in 1:nrow(M)]
    L.ComplexKey = ["" for _ in 1:nrow(L)]
    sqo_mutect = nothing
    if nrow(M) > 0
        sqo_mutect = sqo_main(M,V,"mutect_complex_$(SAMPLE_NAME)",directoryloc)
    end
    if ! isnothing(sqo_mutect)
        sqo_mutect.pass = abs.(log10.(sqo_mutect.VAF ./ Statistics.mean.(sqo_mutect.bestVAF))) .<= 1 .&& [length(i[1])<150 && length(i[2])<150 for i in split.([i[end] for i in split.(sqo_mutect.key," ")],">")] .&& length.(sqo_mutect.calpos_best) .> 1
        perm = sortperm.(sqo_mutect.bestVAF,rev=true)
        chrms = [ i[1] for i in split.(sqo_mutect.key," ")  ]
        m_dict=Dict()
        # for i in 1:nrow(sqo_mutect)
        #     if sqo_mutect[i,"pass"]
        #         s = string(chrms[i]," ", sqo_mutect[i,"calpos_best"][perm[i][1]]," ",sqo_mutect[i,"calref_best"][perm[i][1]],">",sqo_mutect[i,"calalt_best"][perm[i][1]] )
        #         try
        #             l = filter(x->x.key == s,M).Column1[1]
        #             m_dict[ s ] = sqo_mutect[i,"key"]
        #             # push!(M,M[l,:])
        #             # M[l,"key"] = sqo_mutect[i,"key"]
        #             # M[l,"IsComplex"] = true
        #         catch
        #             # l = filter(x->x.key == m_dict[ s ],M).Column1[1]
        #             # push!(M,M[l,:])
        #             # M[end,"key"] = sqo_mutect[i,"key"]
        #             # M[end,"IsComplex"] = true
        #         end    
               
        #     end
        # end
        for i in 1:nrow(sqo_mutect)
            for j in 1:length(perm[i])
                s = string(chrms[i]," ", sqo_mutect[i,"calpos_best"][perm[i][j]]," ",sqo_mutect[i,"calref_best"][perm[i][j]],">",sqo_mutect[i,"calalt_best"][perm[i][j]] )
                m_dict[ s ] = sqo_mutect[i,"key"]
                l = filter(x->x.key == s,M).Column1[1]
                M[l,"IsComplexPart"] = true
                M[l,"ComplexKey"] = sqo_mutect[i,"key"]
            end
            if sqo_mutect[i,"pass"]
                V[filter(x->x.key == sqo_mutect[i,"key"],V).Column1[1],"IsComplex"] = true
                # try
                #     l = filter(x->x.key == s,L).Column1[1]
                #     l_dict[ s ] = sqo_lofreq[i,"key"]
                #     # push!(L,L[l,:])
                #     # L[l,"key"] = sqo_lofreq[i,"key"]
                #     # L[l,"IsComplex"] = true
                # catch
                #     # l = filter(x->x.key == l_dict[ s ],L).Column1[1]
                #     # push!(L,L[l,:])
                #     # L[end,"key"] = sqo_lofreq[i,"key"]
                #     # L[end,"IsComplex"] = true
                # end
            end
        end
    end
    sqo_lofreq = nothing
    if nrow(L) > 0 
        sqo_lofreq = sqo_main(L,V,"lofreq_complex_$(SAMPLE_NAME)",directoryloc)
    end
    if ! isnothing(sqo_lofreq)
        sqo_lofreq.pass = abs.(log10.(sqo_lofreq.VAF ./ Statistics.mean.(sqo_lofreq.bestVAF))) .<= 1 .&& [length(i[1])<150 && length(i[2])<150 for i in split.([i[end] for i in split.(sqo_lofreq.key," ")],">")] .&& length.(sqo_lofreq.calpos_best) .> 1
        perm = sortperm.(sqo_lofreq.bestVAF,rev=true)
        chrms = [ i[1] for i in split.(sqo_lofreq.key," ")  ]
        l_dict=Dict()
        for i in 1:nrow(sqo_lofreq)
            for j in 1:length(perm[i])
                s = string(chrms[i]," ", sqo_lofreq[i,"calpos_best"][perm[i][j]]," ",sqo_lofreq[i,"calref_best"][perm[i][j]],">",sqo_lofreq[i,"calalt_best"][perm[i][j]] )
                l_dict[ s ] = sqo_lofreq[i,"key"]
                l = filter(x->x.key == s,L).Column1[1]
                # println(s," ",l," ",L[l,"key"])
                L[l,"IsComplexPart"] = true
                L[l,"ComplexKey"] = sqo_lofreq[i,"key"]
            end
            if sqo_lofreq[i,"pass"]
                V[filter(x->x.key == sqo_lofreq[i,"key"],V).Column1[1],"IsComplex"] = true
                # try
                #     l = filter(x->x.key == s,L).Column1[1]
                #     l_dict[ s ] = sqo_lofreq[i,"key"]
                #     # push!(L,L[l,:])
                #     # L[l,"key"] = sqo_lofreq[i,"key"]
                #     # L[l,"IsComplex"] = true
                # catch
                #     # l = filter(x->x.key == l_dict[ s ],L).Column1[1]
                #     # push!(L,L[l,:])
                #     # L[end,"key"] = sqo_lofreq[i,"key"]
                #     # L[end,"IsComplex"] = true
                # end
            end
        end
    end

    M.Column1 = 1:nrow(M)
    L.Column1 = 1:nrow(L)
    V.Column1 = 1:nrow(V)

    #### CLEAN ####

    L.sample_key = L.SN_TAG .*" ".* L.key
    M.sample_key = M.SN_TAG .*" ".* M.key
    V.sample_key = V.SN_TAG .*" ".* V.key
    for i in names(L)
        if occursin("lofreq",i)
            s = split(i,"_")
            rename!(L,i=>join(s[1:end-1],"_")*"_Lofreq")
        end
    end
    for i in names(M)
        if occursin("mutect",i)
            s = split(i,"_")
            rename!(M,i=>join(s[1:end-1],"_")*"_Mutect")
        end
    end
    for i in names(V)
        if occursin("vardict",i)
            s = split(i,"_")
            rename!(V,i=>join(s[1:end-1],"_")*"_Vardict")
        end
    end
   
    res = test(M,"FILTER_Mutect")
    M.FILTER_Mutect = res[2]
    M.FP_Filter = res[1]
    M.FPpass = res[3]

    res = test(V,"FILTER_Vardict")
    V.FILTER_Vardict = res[2]
    V.FP_Filter = res[1]
    V.FPpass = res[3]

    res = test(L,"FILTER_Lofreq")
    L.FILTER_Lofreq = res[2];
    L.FP_Filter = res[1];
    L.FPpass = res[3];

    ### RENAME
    rename!(L,:RDF_Lofreq=>:RefFwd);
    rename!(L,:RDR_Lofreq=>:RefRev);
    rename!(L,:ADF_Lofreq=>:AltFwd);
    rename!(L,:ADR_Lofreq=>:AltRev);
    rename!(L,:PON_FISHER=>:pon_pvalue);

    rename!(V,:RDF_Vardict=>:RefFwd);
    rename!(V,:RDR_Vardict=>:RefRev);
    rename!(V,:ADF_Vardict=>:AltFwd);
    rename!(V,:ADR_Vardict=>:AltRev);
    rename!(V,:PON_FISHER=>:pon_pvalue);

    rename!(M,:RDF_Mutect=>:RefFwd);
    rename!(M,:RDR_Mutect=>:RefRev);
    rename!(M,:ADF_Mutect=>:AltFwd);
    rename!(M,:ADR_Mutect=>:AltRev);
    rename!(M,:PON_FISHER=>:pon_pvalue);

    v = filter(x->!(length(x.REF) ==1&&length(x.ALT)==1),V)

    final_df = predict( clean(copy(M),copy(L),copy(V),directoryloc*ponname), modelloc )

    for nam in names(M)
        if ! occursin("utect",nam) && nam != "sample_key" 
            rename!(M,"$(nam)"=>string(nam,"_Mutect_Raw"))
        end
    end

    for nam in names(V)
        if ! occursin("ardict",nam) && nam != "sample_key" 
            rename!(V,"$(nam)"=>string(nam,"_Vardict_Raw")) 
        end
    end

    for nam in names(L)
        if ! occursin("ofreq",nam) && nam != "sample_key" 
            rename!(L,"$(nam)"=>string(nam,"_Lofreq_Raw"))
        end
    end
    
    final = outerjoin(V,L,on=[:sample_key])

    final = outerjoin(final,M,on=[:sample_key]);

    final = leftjoin(final,final_df,on=:sample_key)

    

    final.IsComplex_Mutect_Raw = [ismissing(i) ? false : i for i in final.IsComplex_Mutect_Raw]
    final.IsComplex_Lofreq_Raw = [ismissing(i) ? false : i for i in final.IsComplex_Lofreq_Raw]
    final.IsComplex_Vardict_Raw = [ismissing(i) ? false : i for i in final.IsComplex_Vardict_Raw]
    final.IsComplexPart_Mutect_Raw = [ismissing(i) ? false : i for i in final.IsComplexPart_Mutect_Raw]
    final.IsComplexPart_Lofreq_Raw = [ismissing(i) ? false : i for i in final.IsComplexPart_Lofreq_Raw]
    final.IsComplexPart_Vardict_Raw = [ismissing(i) ? false : i for i in final.IsComplexPart_Vardict_Raw]
    final.ComplexKey_Mutect_Raw = [ismissing(i) ? "" : i for i in final.ComplexKey_Mutect_Raw]
    final.ComplexKey_Lofreq_Raw = [ismissing(i) ? "" : i for i in final.ComplexKey_Lofreq_Raw]
    final.ComplexKey_Vardict_Raw = [ismissing(i) ? "" : i for i in final.ComplexKey_Vardict_Raw]

    pdict = Dict()
    pdict[""] = ""
    pdict[String[]] = ""
    for i in 1:nrow(P)
        pdict[P.key[i]] = ( P.gt_AD_ref[i],P.gt_AD_alt[i])
    end
    
    s = []
    for i in 1:nrow(v)
        kk=filter(x->v.POS[i]-1<x.POS<v.POS[i]+1,P).key
        push!(s,kk[StringDistances.findall(v.key[i], kk, StringDistances.Levenshtein())])
    end
    
    s = [isempty(s[i]) ? [] : (v.key[i],[(k,pdict[k]) for k in s[i]]) for i in 1:nrow(v)]
    filter!(x->!isempty(x),s)
    pindf = DataFrame(key=[i[1] for i in s],PINDEL_MATCH=[i[2] for i in s])
    CSV.write(directoryloc*"output_pindel_complex_$(SAMPLE_NAME).tsv.gz",pindf,delim="\t",compress=true)
    # CSV.write(ARGS[6]*"output_pindel_complex_$(SAMPLE_NAME).tsv.gz",pindf,delim="\t",compress=true)
    final = leftjoin(final,pindf,on=:key)
    # CSV.write("output_$(SAMPLE_NAME).tsv.gz",final,delim="\t",compress=true)


    final.FP_probability = zeros(nrow(final));
    final.pon_FP_pass_XGB = trues(nrow(final))
    final.low_AF_pass_XGB = trues(nrow(final))
    final.long_indel_pass_XGB = trues(nrow(final))
    final.bcbio_pass_XGB = trues(nrow(final))
    final.zscore_pass_XGB = trues(nrow(final))
    final.all_fp_pass_XGB = trues(nrow(final))
    final.di_tri_vard_pass_XGB=trues(nrow(final));
    final.long100_indel_pass_XGB = trues(nrow(final))
    final.PASS_BY_1 = falses(nrow(final))
    final.PASS_BY_2 = falses(nrow(final))
    final.PASS_BY_3 = falses(nrow(final))
    for i in 1:nrow(final)
        # if ! isnan(x.int_pon_pvalue_lofreq[i]) && x.int_pon_pvalue_lofreq[i] > -5.673835246901769 #log10(2.119164905e-6) # -5.673835246901769
        #        x.FP_probability[i] -= 1.11
        #    elseif ! isnan(x.int_pon_pvalue_mutect[i]) && x.int_pon_pvalue_mutect[i] > -5.673835246901769 # log10(2.119164905e-6)
        #        x.FP_probability[i] -= 1.12
        #    elseif ! isnan(x.int_pon_pvalue_vardict[i]) && x.int_pon_pvalue_vardict[i] > -5.673835246901769 # log10(2.119164905e-6)
        #        x.FP_probability[i] -= 1.14
        #    end
       passbyv = final.FILTER_Lofreq_PASS_XGB[i] + final.FILTER_Vardict_PASS_XGB[i] + final.FILTER_Mutect_PASS_XGB[i]
       if passbyv > 2.5
        final.PASS_BY_3[i] = true
       end
       if passbyv > 1.5
        final.PASS_BY_2[i] = true
       end
       if passbyv > 0.5
        final.PASS_BY_1[i] = true
       end
   
       if ! isnan(final.pon_pvalue_Lofreq_XGB[i]) && final.pon_pvalue_Lofreq_XGB[i] > log10(2.119164905e-6) # -5.67383
            final.FP_probability[i] -= 2.11
            final.pon_FP_pass_XGB[i] = false
       elseif ! isnan(final.pon_pvalue_Mutect_XGB[i]) && final.pon_pvalue_Mutect_XGB[i] > log10(2.119164905e-6)
            final.FP_probability[i] -= 2.12
            final.pon_FP_pass_XGB[i] = false
       elseif ! isnan(final.pon_pvalue_Vardict_XGB[i]) && final.pon_pvalue_Vardict_XGB[i] > log10(2.119164905e-6)
            final.FP_probability[i] -= 2.14
            final.pon_FP_pass_XGB[i] = false
       end
   
       if final.AF_XGB[i] < 0.001
            final.FP_probability[i] -= 4.0
           final.low_AF_pass_XGB[i] = false
       end

       if final.FPpass_XGB[i] == false
        final.FP_probability[i] -= 8.0
       end

       if final[i,:ref_len_XGB] > 100 || final[i,:alt_len_XGB] > 100
            final[i,:FP_probability] -= 16.5
            final.long100_indel_pass_XGB[i] = false
       end
   
       if (final[i,:ref_len_XGB] > 20 || final[i,:alt_len_XGB] > 20) && final[i,:FILTER_Lofreq_Not_Detected_XGB] > 0.5 && final[i,:FILTER_Mutect_Not_Detected_XGB] > 0.5 && ismissing(final[i,:PINDEL_MATCH])
            #  final[i,:FP_probability] -= 16.0
            final.long_indel_pass_XGB[i] = false
       end

       if final[i,:ref_len_XGB] == final[i,:alt_len_XGB] && final[i,:alt_len_XGB]>1 && final[i,:FILTER_Lofreq_Not_Detected_XGB] > 0.5 && final[i,:FILTER_Mutect_Not_Detected_XGB] > 0.5 && ismissing(final[i,:PINDEL_MATCH])
            final[i,:FP_probability] -= 1024.0
            final.di_tri_vard_pass_XGB[i] = false
        end
       
       if final[i,:FILTER_Vardict_BCBIO_XGB] > 0.5
        final[i,:FP_probability] -= 32.0
           final.bcbio_pass_XGB[i] = false
       end
       
       if final.pon_af_zscore_XGB[i] < 1.96
           # df2[i,"Real"] = false
           final[i,:FP_probability] -= 512.0
           final.zscore_pass_XGB[i] = false
       end
    #    if final[i,:FP_probability] < -0.01
    #     final.all_fp_pass_XGB[i] = false
    # end
   end

    #    for i in ["pon_FP_pass","low_AF_pass","long_indel_pass","bcbio_pass","zscore_pass","all_fp_pass","di_tri_vard_pass"]
    #         rename!(x,i=>string(i,"_XGB"))
    #     end

    select!(final,Not(["key_Vardict_Raw","key_Lofreq_Raw","key_Mutect_Raw","Indiv_Vardict_Raw","Indiv_Mutect_Raw",
    "gt_AD_ref_Vardict_XGB","gt_AD_ref_Lofreq_XGB","gt_AD_ref_Mutect_XGB",
    # "FP_probability",
    "QUAL_Vardict_XGB",
    "PMEAN_Vardict_XGB","ODDRATIO_Vardict_XGB","SN_Vardict_XGB","SHIFT3_Vardict_XGB","NM_Vardict_XGB",
    "HICNT_Vardict_XGB","HICOV_Vardict_XGB","QUAL_Lofreq_XGB","GERMQ_Mutect_XGB",
    "MPOS_Mutect_XGB","ROQ_Mutect_XGB","TLOD_Mutect_XGB"
     ]))
    
    nms=names(final)
    
    # nmsr = replace.(nms,"_Raw"=>"")
    nmsr = replace.(nms,"_Vardict_Raw"=>"")
    nmsr = replace.(nmsr,"_Lofreq_Raw"=>"")
    nmsr = replace.(nmsr,"_Mutect_Raw"=>"")
    
    # println(names(final)[occursin.("omplex",names(final))])

    combinerraw!(final,"IsComplex")
    for (k,v) in StatsBase.countmap(nmsr)
        if v>1 && ! occursin("PON",k) && ! occursin("pon",k) && ! occursin("gt_",k) && ! occursin("AltRev",k) && ! occursin("RefRev",k) &&
            ! occursin("AltFwd",k) && ! occursin("RefFwd",k) && ! occursin("DP",k) && "IsComplex" != k
            # println(k) #," ",v)
            combinerraw!(final,k)
        end
    end
    # println(names(final)[occursin.("omplex",names(final))])
    # allowmissing!(final)
    # for j in 1:ncol(final), i in 1:nrow(final)
    #     if ! ismissing(final[i,j]) && (final[i,j] == NaN32 || final[i,j] == NaN)
    #         final[i,j] = missing
    #     end
    # end

    mcomp=CSV.read(directoryloc*"output_mutect_complex_$(SAMPLE_NAME).tsv.gz",DataFrame,delim="\t",downcast=false,stringtype=String, missingstring=["","NA", "NAN", "NULL"],ntasks=1,truestrings=["true", "True", "TRUE", "1"],falsestrings=["false", "False", "FALSE", "0"]);
    lcomp=CSV.read(directoryloc*"output_lofreq_complex_$(SAMPLE_NAME).tsv.gz",DataFrame,delim="\t",downcast=false,stringtype=String, missingstring=["","NA", "NAN", "NULL"],ntasks=1,truestrings=["true", "True", "TRUE", "1"],falsestrings=["false", "False", "FALSE", "0"]);
    rename!(mcomp,"calpos_best"=>"complex_mutect_pos")
    rename!(mcomp,"calref_best"=>"complex_mutect_ref")
    rename!(mcomp,"calalt_best"=>"complex_mutect_alt")
    rename!(mcomp,"bestVAF"=>"complex_mutect_vaf")
    rename!(lcomp,"calpos_best"=>"complex_lofreq_pos")
    rename!(lcomp,"calref_best"=>"complex_lofreq_ref")
    rename!(lcomp,"calalt_best"=>"complex_lofreq_alt")
    rename!(lcomp,"bestVAF"=>"complex_lofreq_vaf")
    select!(mcomp,Not("VAF"))
    select!(lcomp,Not("VAF"))
    final = leftjoin(final,lcomp,on=:key)
    final = leftjoin(final,mcomp,on=:key)

    final.PossibleComplex = falses(nrow(final))
    for i in 1:nrow(final)
        if ! ismissing(final.complex_mutect_vaf[i]) || ! ismissing(final.complex_lofreq_vaf[i]) 
            final.PossibleComplex[i] = true
        end
    end

    for i in 1:nrow(final) 
        if ! final.long_indel_pass_XGB[i] 
            if final.PossibleComplex[i]
                final.long_indel_pass_XGB[i] = true
            else
                final[i,:FP_probability] -= 16.0
            end
        end
      
        if final[i,:FP_probability] < -0.1
            final.all_fp_pass_XGB[i] = false
        end
    end

    finalcolumnorder = [
        "CHROM","POS","ID","REF","ALT","key","sample_key","subject","Gene","OLD_MULTIALLELIC","SN_TAG","QUAL_Mutect","FILTER_Mutect","AS_FilterStatus_Mutect_Raw",
        "AS_SB_TABLE_Mutect_Raw","AS_UNIQ_ALT_READ_COUNT_Mutect_Raw","CONTQ_Mutect_Raw","DP_Mutect_Raw","ECNT_Mutect_Raw","GERMQ_Mutect_Raw",
        "MBQ_Mutect_Raw","MFRL_Mutect_Raw","MMQ_Mutect_Raw","MPOS_Mutect_Raw","NALOD_Mutect_Raw","NCount_Mutect_Raw","NLOD_Mutect_Raw","OCM_Mutect_Raw",
        "PON_Mutect_Raw","POPAF_Mutect_Raw","ROQ_Mutect_Raw","RPA_Mutect_Raw","RU_Mutect_Raw","SEQQ_Mutect_Raw","STR_Mutect_Raw","STRANDQ_Mutect_Raw",
        "STRQ_Mutect_Raw","TLOD_Mutect_Raw","PON_2AT2_percent_Mutect_Raw","PON_NAT2_percent_Mutect_Raw","PON_MAX_VAF_Mutect_Raw","gt_AD_ref_Mutect",
        "gt_AD_alt_Mutect","gt_AF_Mutect","gt_DP_Mutect","gt_F1R2_Mutect","gt_F2R1_Mutect","gt_GQ_Mutect","gt_GT_Mutect","gt_PGT_Mutect",
        "gt_PID_Mutect","gt_PL_Mutect","gt_PS_Mutect","RefFwd_Mutect_Raw","RefRev_Mutect_Raw","AltFwd_Mutect_Raw","AltRev_Mutect_Raw","gt_GT_alleles_Mutect",
        "pass_strand_bias_Mutect","QUAL_Lofreq","FILTER_Lofreq","DP_Lofreq_Raw","gt_AF_Lofreq","SB_Lofreq_Raw","RefFwd_Lofreq_Raw","RefRev_Lofreq_Raw",
        "AltFwd_Lofreq_Raw","AltRev_Lofreq_Raw","INDEL_Lofreq_Raw","CONSVAR_Lofreq_Raw","HRUN_Lofreq_Raw","PON_2AT2_percent_Lofreq_Raw","PON_NAT2_percent_Lofreq_Raw",
        "PON_MAX_VAF_Lofreq_Raw","gt_AD_alt_Lofreq","gt_AD_ref_Lofreq","pass_strand_bias_Lofreq","QUAL_Vardict","FILTER_Vardict","TYPE_Vardict_Raw","DP_Vardict_Raw",
        "END_Vardict_Raw","VD_Vardict_Raw","AF_Vardict_Raw","BIAS_Vardict_Raw","REFBIAS_Vardict_Raw","VARBIAS_Vardict_Raw","PMEAN_Vardict_Raw","PSTD_Vardict_Raw",
        "ReadQual_Vardict_Raw","QSTD_Vardict_Raw","SBF_Vardict_Raw","ODDRATIO_Vardict_Raw","MQ_Vardict_Raw","SN_Vardict_Raw","HIAF_Vardict_Raw","ADJAF_Vardict_Raw",
        "SHIFT3_Vardict_Raw","MSI_Vardict_Raw","MSILEN_Vardict_Raw","NM_Vardict_Raw","LSEQ_Vardict_Raw","RSEQ_Vardict_Raw","GDAMP_Vardict_Raw","TLAMP_Vardict_Raw",
        "NCAMP_Vardict_Raw","AMPFLAG_Vardict_Raw","HICNT_Vardict_Raw","HICOV_Vardict_Raw","SPLITREAD_Vardict_Raw","SPANPAIR_Vardict_Raw","SVTYPE_Vardict_Raw",
        "SVLEN_Vardict_Raw","DUPRATE_Vardict_Raw","PON_2AT2_percent_Vardict_Raw","PON_NAT2_percent_Vardict_Raw","PON_MAX_VAF_Vardict_Raw","gt_GT_Vardict",
        "gt_DP_Vardict","gt_VD_Vardict","gt_AD_ref_Vardict","gt_AD_alt_Vardict","gt_AF_Vardict","RefFwd_Vardict_Raw","RefRev_Vardict_Raw",
        "AltFwd_Vardict_Raw","AltRev_Vardict_Raw","gt_GT_alleles_Vardict","pass_strand_bias_Vardict","Allele_VEP","Consequence_VEP",
        "IMPACT_VEP","SYMBOL_VEP","Gene_VEP","Feature_type_VEP","Feature_VEP","BIOTYPE_VEP","EXON_VEP","INTRON_VEP","cDNA_position_VEP",
        "CDS_position_VEP","Protein_position_VEP","Amino_acids_VEP","Codons_VEP","Existing_variation_VEP","DISTANCE_VEP",
        "STRAND_VEP","FLAGS_VEP","VARIANT_CLASS_VEP","SYMBOL_SOURCE_VEP","HGNC_ID_VEP","CANONICAL_VEP","MANE_VEP",
        "TSL_VEP","APPRIS_VEP","CCDS_VEP","ENSP_VEP","SWISSPROT_VEP","TREMBL_VEP","UNIPARC_VEP","REFSEQ_MATCH_VEP","SOURCE_VEP",
        "REFSEQ_OFFSET_VEP","GIVEN_REF_VEP","USED_REF_VEP","BAM_EDIT_VEP","GENE_PHENO_VEP","SIFT_VEP","PolyPhen_VEP",
        "DOMAINS_VEP","miRNA_VEP","HGVS_OFFSET_VEP","CLIN_SIG_VEP","SOMATIC_VEP","PHENO_VEP","PUBMED_VEP","VAR_SYNONYMS_VEP",
        "MOTIF_NAME_VEP","MOTIF_POS_VEP","HIGH_INF_POS_VEP","MOTIF_SCORE_CHANGE_VEP","TRANSCRIPTION_FACTORS_VEP",
        "FrameshiftSequence_VEP","WildtypeProtein_VEP","AF_VEP","AFR_AF_VEP","AMR_AF_VEP","EAS_AF_VEP","EUR_AF_VEP","SAS_AF_VEP",
        "AA_AF_VEP","EA_AF_VEP","gnomAD_AF_VEP","gnomAD_AFR_AF_VEP","gnomAD_AMR_AF_VEP","gnomAD_ASJ_AF_VEP","gnomAD_EAS_AF_VEP",
        "gnomAD_FIN_AF_VEP","gnomAD_NFE_AF_VEP","gnomAD_OTH_AF_VEP","gnomAD_SAS_AF_VEP","MAX_AF_VEP","MAX_AF_POPS_VEP","gnomADe_VEP","gnomADe_AF_VEP",
        "gnomADe_AF_AFR_VEP","gnomADe_AF_AMR_VEP","gnomADe_AF_ASJ_VEP","gnomADe_AF_EAS_VEP","gnomADe_AF_FIN_VEP","gnomADe_AF_NFE_VEP","gnomADe_AF_OTH_VEP",
        "gnomADe_AF_SAS_VEP","gnomADg_VEP","gnomADg_AF_VEP","gnomADg_AF_ami_VEP","gnomADg_AF_oth_VEP","gnomADg_AF_afr_VEP","gnomADg_AF_sas_VEP","gnomADg_AF_asj_VEP",
        "gnomADg_AF_fin_VEP","gnomADg_AF_amr_VEP","gnomADg_AF_nfe_VEP","gnomADg_AF_eas_VEP","gnomADe","gnomADe_AF","gnomADe_AF_AFR","gnomADe_AF_AMR","gnomADe_AF_ASJ",
        "gnomADe_AF_EAS","gnomADe_AF_FIN","gnomADe_AF_NFE","gnomADe_AF_OTH","gnomADe_AF_SAS","gnomADg","gnomADg_AF","gnomADg_AF_ami","gnomADg_AF_oth","gnomADg_AF_afr",
        "gnomADg_AF_sas","gnomADg_AF_asj","gnomADg_AF_fin","gnomADg_AF_amr","gnomADg_AF_nfe","gnomADg_AF_eas","clinvar_VEP","clinvar_CLINSIGN_VEP","clinvar_PHENOTYPE_VEP",
        "clinvar_SCORE_VEP","clinvar_RCVACC_VEP","clinvar_TESTEDINGTR_VEP","clinvar_PHENOTYPELIST_VEP","clinvar_NUMSUBMIT_VEP","clinvar_GUIDELINES_VEP",
        "clinvar","clinvar_CLINSIGN","clinvar_PHENOTYPE","clinvar_SCORE","clinvar_RCVACC","clinvar_TESTEDINGTR","clinvar_PHENOTYPELIST","clinvar_NUMSUBMIT",
        "clinvar_GUIDELINES","MFRL_Mutect_2_XGB","MMQ_Mutect_1_XGB","MMQ_Mutect_2_XGB","FPpass_XGB","PONvafmean_XGB","PONvafstd_XGB","pon_af_zscore_XGB",
        "FILTER_Mutect_Not_Detected_XGB","FILTER_Mutect_clustered_events_XGB","FILTER_Mutect_haplotype_XGB","FILTER_Mutect_orientation_XGB","FILTER_Mutect_position_XGB",
        "FILTER_Mutect_strand_bias_XGB","FILTER_Mutect_PASS_XGB","FILTER_Mutect_weak_evidence_XGB","FILTER_Lofreq_Not_Detected_XGB","FILTER_Lofreq_PASS_XGB",
        "FILTER_Vardict_PASS_XGB","FILTER_Vardict_NM5_25_XGB","FILTER_Vardict_p8_XGB","FILTER_Vardict_pSTD_XGB","FILTER_Vardict_BCBIO_XGB",
        "FILTER_Vardict_Not_Detected_XGB","FP_Filter_PASS_XGB","FP_Filter_RLD25_XGB","FP_Filter_DETP20_XGB","FP_Filter_MVC4_XGB","FP_Filter_SB1_XGB",
        "FP_Filter_MMQSD50_XGB","FP_Filter_NRC_XGB","FP_Filter_PB10_XGB","FP_Filter_MMQS100_XGB","FPpass","FP_Filter","dustscore_XGB",
        "dustscore3_XGB","dustscore5_XGB","dustscore10_XGB","passhomopolymerfilter_XGB","maxgnomADAF_XGB","PONRefDepth_XGB","PONAltDepth_XGB","pon_pvalue_Mutect_Raw",
        "pon_pvalue_Lofreq_Raw","pon_pvalue_Vardict_Raw","pon_pvalue_Mutect_XGB","pon_pvalue_Lofreq_XGB","pon_pvalue_Vardict_XGB","clinvarSCOREVEP_XGB",
        "sourcetotalsloci_XGB","sourcetotalsp_XGB","sourcetotalsc_XGB","CosmicCount_XGB","hemecosmiccount_XGB","myeloidcosmiccount_XGB","predicted_probability_XGB",
        "AF_XGB","PINDEL_MATCH","pon_FP_pass_XGB","low_AF_pass_XGB","long_indel_pass_XGB","bcbio_pass_XGB","zscore_pass_XGB","all_fp_pass_XGB","di_tri_vard_pass_XGB",
        "long100_indel_pass_XGB","max_gnomAD_AF","max_gnomADe_AF_VEP","max_gnomADg_AF_VEP","max_gnomAD_AF_VEP","pass_homopolymer_filter","case_NXXX","case_XNXX",
        "case_XXNX","case_XXXN","case_NNXX","case_XNNX","case_XXNN","context_3","context_10","dust_score_5","dust_score_3","dust_score_10","dust_score",
        "AAchange.x","gene_loci_p","gene_loci_c","gene_loci_vep","gene_aachange","gene_cDNAchange","n.loci.vep","source.totals.loci","n.HGVSp.x","source.totals.p","n.HGVSc.x",
        "source.totals.c","n.HGVSp.y","n.HGVSc.y","COSMIC_ID","CosmicCount","heme_cosmic_count","myeloid_cosmic_count","oncoKB","isOncogenic","isTSG","isTruncatingHotSpot","ch_my_pd",
        "ch_pd","ch_pd2","VariantClass","AAchange.y","WHY_CH","PASS_BY_1","PASS_BY_2","PASS_BY_3","PossibleComplex","IsComplexPart","ComplexKey",
        "complex_lofreq_pos","complex_lofreq_ref","complex_lofreq_alt",
        "complex_lofreq_vaf","complex_mutect_pos","complex_mutect_ref","complex_mutect_alt","complex_mutect_vaf"
    ]

    


    CSV.write(
        directoryloc*"output_$(SAMPLE_NAME).tsv.gz",
        select(final,finalcolumnorder),
        delim="\t",
        compress=true
        )
    return nothing
end

function real_main()::Int
    main(ARGS[1],ARGS[2],ARGS[3],ARGS[4],ARGS[5],ARGS[6],"/opt/bin/XGB_FINAL_MODEL.jlso")
    # main(
    #     "lofreq.218281_T_S209.final.annotated.tsv",
    #     "mutect.218281_T_S209.final.annotated.tsv",
    #     "vardict.218281_T_S209.final.annotated.tsv",
    #     "pindel_full.218281_T_S209.vcf.gz",
    #     "218281_T_S209.pon.total.counts.vcf.gz",
    #     "",
    #     "XGB_FINAL_MODEL.jlso"
    # )
    return 0
end
