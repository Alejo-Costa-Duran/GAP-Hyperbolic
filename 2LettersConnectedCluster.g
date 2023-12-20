# 2LettersConnectedCluster.g
# This program computes the PBC hopping matrices of clusters of the {14,7} lattice
# of size n=16,...,182, belonging to the 2 letter words generation

# this version first ensures that the length-1 words appear as distinct elements in the tranversal

# assumes that LINS has been loaded and NSGbolza.g run, storing the normal subgroups in the list L

# read all 56 nontrivial length-2 elements (enumerated using Matlab)
# and find which ones appear in the transversal
# maximal index that LINS has computed
nmax:=23;;

# index for which we want to compute the hopping matrices and output the transversal
n:=23;;


Read("./2LetterWords.g");


# initialize empty file: transversal data
filenameT:="./PBCData/p2anisT_n17.m";;
PrintTo(filenameT);;
# initialize empty file: hopping matrix
filenameH:="./PBCData/p2anisH_n17.m";;
PrintTo(filenameH);;
# initialize empty file: translation matrices
# (only root of filename here, subgroup label will be appended later)
rootfnameU:="./PBCData/p2anisU_n17_k";;

belown:=Sum(nsg{[1..n-1]});;
# loop through all subgroups of index n
kanis:=0;;
p2anis:=[];;
for k in [1..100] do
	# construct transversal T_n^{(k)}, k = 1,...,n
	T:=RightTransversal(G,Grp(List_NSG[belown+k]));;
	# canonical position of cosets
	p2:=[];;
    	p2[1]:=PositionCanonical(T,One(G));;
	# 14 length-1 words
	for m in [1..14] do
        	p2[1+m]:=PositionCanonical(T,gen[m]);;
    od;
    for m in [2..15] do
        Print(p2[m],"\t");
        od;
        Print("\n");
	# check that identity is the first element	
	if p2[1] <> 1 then
        	Print("First coset not identity for k=",k,"!\n");
        	break;
	# check that length-1 words appear as distinct elements in the transversal
    	elif IsDuplicateFree(p2{[1..14]}) then
            Print("Is Duplicate free");
        	for m in [1..182] do # 182 length-2 words
			p2[15+m]:=PositionCanonical(T,p2wds[m]);
		od;
        	# search for clusters: array p2 contains {1,2,...,n}
        	if IsSubsetSet(p2,[1..n]) then
            		Print("found! L[",belown+k,"].Group\n");
            		kanis:=kanis+1;;
            		p2anis[kanis]:=belown+k;;
            		# create list of coset representatives g_i in canonical order
            		g:=[];;
			glabel:=[];; # group element identifier, for Matlab purposes
            		# do length-2 words first, some of which might overlap with
            		# length-1 words in the transversal
            		for m in [1..182] do
                		# if two length-2 words belong to the same coset, the
                		# last will be chosen (arbitrary choice)
                		g[p2[15+m]]:=p2wds[m];;
				glabel[p2[15+m]]:=14+m;; # entry 9 and above means length-2 word
            		od;
            		# length-1 words done afterwards, to ensure they appear in the
            		# choice of transversal
            		for m in [1..14] do
                		g[p2[1+m]]:=gen[m];;
				glabel[p2[1+m]]:=m;; # entry 1-8 means length-1 word 
            		od;
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
            		p2test:=[];;
            		for i in [1..n] do
                		p2test[i]:=PositionCanonical(T,g[i]);;
            		od;
			# array g[i] should contain elements of the transversal
			# in canonical order
            		if p2test <> [1..n] then
                		Print("Error #3 in constructing transversal!\n");
                		break;
            		fi;
			# output transversal data to file
			AppendTo(filenameT,"SgLabel(",kanis,")=",p2anis[kanis],";\n");			
			for i in [1..n] do
				AppendTo(filenameT,"glabel(",i,",",kanis,")=",glabel[i],";\n");
			od;
            		# initialize hopping matrix to zero
            		H:=DiagonalMat([1..n]);;
            		for i in [1..n] do
                		H[i][i]:=0;;
            		od;
            		for i in [1..n] do # pick a site
                		for a in [1..14] do # find its 8 nearest neighbors
                    			j:=PositionCanonical(T,g[i]*gen[a]);
                    			if i <> j then
                        			H[i][j]:=-1;;
                    			fi;
                		od;
            		od;
            		# output hopping matrices in Matlab-ready format
            		AppendTo(filenameH,"SgLabel(",kanis,")=",p2anis[kanis],";\n");
            		AppendTo(filenameH,"H(:,:,",kanis,")=[");;
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
			# construct translation U matrices
			filenameU:=Concatenation(rootfnameU,PrintString(kanis),".m");;
			PrintTo(filenameU);;
			for kk in [1..n] do # U_{ij}(g_k) for each k=1,...,n
				# initialize U matrix to zero
				U:=DiagonalMat([1..n]);;
				for i in [1..n] do
					U[i][i]:=0;;
				od;
				for j in [1..n] do
					i:=PositionCanonical(T,g[kk]*g[j]);;
					U[i][j]:=1;;
				od;
				# output U matrix in Matlab-ready format
				AppendTo(filenameU,"U(:,:,",kk,")=[");;
                		for i in [1..n-1] do # write first n-1 rows of U matrix
                        		for j in [1..n] do
                                		AppendTo(filenameU," ",U[i][j]);;
                        		od;
                        		AppendTo(filenameU,";\n");;
                		od;
                		for j in [1..n] do # write last row of U matrix
                        		AppendTo(filenameU," ",U[n][j]);;
                		od;
                		AppendTo(filenameU,"];\n");;
			od;
        	fi;
        fi;

od;
Print("Found ",kanis," anisotropic clusters with ",n," sites\n");


	
