# PBCp1klein.g

# This program constructs connected PBC clusters of the {14,7} lattice
# of size n=2,...,15, belonging to the p=1 generation 

# To the difference of PBCp1bolza.g and PBCp2bolza.g for the {8,8} lattice,
# this program outputs the transversal, hopping matrix, translation matrices, AND the coset table.
# It also outputs the type (abelian/nonabelian) of cover.

# Here, we use the original Klein side-pairings (NOT opposite-side identification)

# assumes that LINS has been loaded and NSGklein.g run, storing the normal subgroups in the list L

# maximal index that LINS has computed
nmax:=4;;

# index for which we want to compute everything (cannot exceed 15)
n:=4;;

# initialize empty file: transversal data
filenameT:=Concatenation("~/Documents/HyperbolicLatticesData/TransversalInfo",String(n),"Clusters.py");;
PrintTo(filenameT);;
# initialize empty file: hopping matrix
filenameH:="~/Documents/HyperbolicLatticesData/p1kleinH_klein_n8.py";;
PrintTo(filenameH);;
# initialize empty file: translation matrices
# (only root of filename here, subgroup label will be appended later)
rootfnameU:="~/Documents/HyperbolicLatticesData/p1kleinU_klein_n8_k";;
# initialize empty file: coset table of factor group G/G_{PBC}
filenameCT:=Concatenation("~/Documents/HyperbolicLatticesData/CosetTable",String(n),"Clusters.py");;
PrintTo(filenameCT);;
# initialize empty file: cover type (abelian/nonabelian)
filenameC:="~/Documents/HyperbolicLatticesData/p1kleinC_klein_n8.py";;
PrintTo(filenameC);;

belown:=Sum(nsg{[1..n-1]});;
# loop through all subgroups of index n
kklein:=0;;
p1klein:=[];;
for k in [1..nsg[n]] do
	# construct transversal T_n^{(k)}, k = 1,...,n
	T:=RightTransversal(G,Grp(List(List_NSG)[belown+k]));;
	# canonical position of cosets
	p1:=[];;
    p1[1]:=PositionCanonical(T,One(G));;
	# 14 length-1 words
	for m in [1..14] do
        	p1[1+m]:=PositionCanonical(T,gen[m]);;
	od;
	Print(p1);
	# check that identity is the first element	

        Print("found! L[",belown+k,"].Group\n");
        kklein:=kklein+1;;
        p1klein[kklein]:=belown+k;;
        # create list of coset representatives g_i in canonical order
        g:=[];;
		glabel:=[];; # group element identifier, for Matlab purposes
        for m in [1..14] do
		# if two length-1 words belong to the same coset, the
		# last will be chosen (arbitrary choice)
            g[p1[1+m]]:=gen[m];;
			glabel[p1[1+m]]:=m;; # entry 1-14 means length-1 word 
        od;
		Print(g);
        # finally, ensure identity element hasn't been rewritten
        g[1]:=One(G);;
		glabel[1]:=0;; # 0 means identity element
            	# ensure all elements of the tranversal are distinct
        if IsDuplicateFree(g) = false then
          	Print("Error #1 in constructing transversal!\n");
           	break;
        fi;
		if IsDuplicateFree(glabel) = false then
			Print("Error #2 in constructing transversal!\n");
			break;
		fi;
        p1test:=[];;
        for i in [1..n] do
           	p1test[i]:=PositionCanonical(T,g[i]);;
        od;
		# array g[i] should contain elements of the transversal
		# in canonical order
        if p1test <> [1..n] then
           	Print("Error #3 in constructing transversal!\n");
           	break;
        fi;
		# output cover type (abelian/nonabelian)
		if IsAbelian(FactorGroup(G,Grp(List(List_NSG)[p1klein[kklein]]))) then # abelian
			AppendTo(filenameC,"% L[",p1klein[kklein],"].Group:\n");
			AppendTo(filenameC,"abelian(",kklein,")=1;\n");
		else # nonabelian
			AppendTo(filenameC,"% L[",p1klein[kklein],"].Group:\n"); 
            AppendTo(filenameC,"abelian(",kklein,")=0;\n"); 	
		fi;
		# output transversal data to file
		AppendTo(filenameT,"SubGroupLabel[",kklein,"]=",p1klein[kklein],"\n");			
		for i in [1..n] do
			AppendTo(filenameT,"glabel[",i,",",kklein,"]=",glabel[i],"\n");
		od;
		# ga(alpha) = coset to which generator \gamma_\alpha belongs (\alpha=1,...,14)	
		AppendTo(filenameT,"ga[",kklein,"]=[");
		for m in [1..13] do
			AppendTo(filenameT,p1[1+m]-1,",");
		od;
		AppendTo(filenameT,p1[1+14]-1,"]\n");
		# compute coset table
		CT:=DiagonalMat([1..n]);;
                for i in [1..n] do
                	CT[i][i]:=0;;
                od;
                for i in [1..n] do
                	for j in [1..n] do
                        	CT[i][j]:=PositionCanonical(T,g[i]*g[j]);;
                        od;
                od;
                # output coset table in Matlab format
                AppendTo(filenameCT,"# L[",p1klein[kklein],"].Group\n");
                AppendTo(filenameCT,"CT[:,:,",kklein,"]=[");;
                for i in [1..n-1] do # write first n-1 rows of coset table
                	AppendTo(filenameCT,"[",CT[i][1]-1);;
					for j in [2..n-1] do
                    	AppendTo(filenameCT,",",CT[i][j]-1);;
                    od;
                    AppendTo(filenameCT,",",CT[i][n]-1,"],");;
                od;
                AppendTo(filenameCT,"[",CT[n][1]-1);;
				for j in [2..n-1] do
                   	AppendTo(filenameCT,",",CT[n][j]-1);;
                od;
                AppendTo(filenameCT,",",CT[n][n]-1,"]");;
                AppendTo(filenameCT,"]\n");;
            	# initialize hopping matrix to zero
            	H:=DiagonalMat([1..n]);;
            	for i in [1..n] do
                	H[i][i]:=0;;
            	od;
            	for i in [1..n] do # pick a site
                	for a in [1..n] do # find its 30 nearest neighbors
                    		j:=PositionCanonical(T,g[i]*gen[a]);
                        	H[i][j]:=H[i][j]+1;;
                	od;
            	od;
            	# output hopping matrices in Matlab-ready format
            	AppendTo(filenameH,"SgLabel(",kklein,")=",p1klein[kklein],";\n");
            	AppendTo(filenameH,"H(:,:,",kklein,")=[");;
            	for i in [1..n-1] do # write first n-1 rows of H matrix
                	for j in [1..n] do
                    		AppendTo(filenameH," ",H[i][j]);;
                	od;
                	AppendTo(filenameH,";\n");;
            	od;
            	for j in [1..n] do # write last row of H matrix
                	AppendTo(filenameH," ",H[n][j]);;
            	od;
            	AppendTo(filenameH,"];\n");;

od;
Print("Found ",kklein," clusters with ",n," sites\n");