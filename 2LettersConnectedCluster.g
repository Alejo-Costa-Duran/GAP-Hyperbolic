# 2LettersConnectedCluster.g
# This program computes the PBC hopping matrices of clusters of the {14,7} lattice
# of size n=16,...,182, belonging to the 2 letter words generation

# this version first ensures that the length-1 words appear as distinct elements in the tranversal

# assumes that LINS has been loaded and NSGbolza.g run, storing the normal subgroups in the list L

# read all 56 nontrivial length-2 elements (enumerated using Matlab)
# and find which ones appear in the transversal
# maximal index that LINS has computed
nmax:=4;;

# index for which we want to compute the hopping matrices and output the transversal
n:=4;;


Read("./2LetterWords.g");


# initialize empty file: transversal data
filenameT:=Concatenation("../PBCData/TransversalInfo,",String(n),"Clusters.py");;
PrintTo(filenameT);;
# initialize empty file: hopping matrix
filenameCT:=Concatenation("../PBCData/CosetTable",String(n),"Clusters.py");;


AppendTo(filenameT,"import numpy as np \nSubgroupLabel = np.zeros(",nsg[n]+1,") \n", "glabel = np.zeros((",n+1,",",nsg[n]+1,")) \n", "ga = [[] for n in range(",nsg[n],")] \n");

AppendTo(filenameCT,"import numpy as np\nCT = np.zeros((",n,",",n,",",nsg[n]+1,"))\n");
belown:=Sum(nsg{[1..n-1]});;
# loop through all subgroups of index n
kanis:=0;;
p2anis:=[];;
for k in [1..nsg[n]] do
	# construct transversal T_n^{(k)}, k = 1,...,n
	T:=RightTransversal(G,Grp(List_NSG[belown+k]));;
	# canonical position of cosets
	p2:=[];;
    	p2[1]:=PositionCanonical(T,One(G));;
	# 14 length-1 words
	for m in [1..14] do
        	p2[1+m]:=PositionCanonical(T,gen[m]);;
    od;
	for m in [1..182] do # 182 length-2 words
			p2[15+m]:=PositionCanonical(T,p2wds[m]);

		od;
	# check that identity is the first element	
	if p2[1] <> 1 then
        	Print("First coset not identity for k=",k,"!\n");
        	break;
	fi;
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
		AppendTo(filenameT,"SubGroupLabel[",kklein,"]=",p1klein[kklein],"\n");			
		for i in [1..n] do
			AppendTo(filenameT,"glabel[",i,",",kklein,"]=",glabel[i],"\n");
		od;
		AppendTo(filenameT,"ga[",kanis,"]=[");
		for m in [1..13] do
			AppendTo(filenameT,p2[1+m]-1,",");
		od;
			AppendTo(filenameT,p2[1+14]-1,"]\n");
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
                AppendTo(filenameCT,"# L[",p2anis[kanis],"].Group\n");
                AppendTo(filenameCT,"CT[:,:,",kanis,"]=[");;
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
 od;
Print("Found ",kanis," anisotropic clusters with ",n," sites\n");