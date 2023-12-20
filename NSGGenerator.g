# NSGGenerator.g
# This program runs LINS to compute the first n_group normal subgroups of the Fuchsian group of the {14,7} lattice
# up to index n_max

# here, use Klein identification (NOT opposite-side)
# presentation with 7 generators + 2 relators

# load low-index normal subgroup package by Friedrich Rober
# (GAP implementation of the algorithm developed by D. Firth and D. Holt for Magma)
Print("Loading LINS package... ");
LoadPackage("LINS");;
Print("done\n");

Print("Defining fundamental group... ");
# create free group F on g1, g2, g3, g4, g5, g6, g7
F:=FreeGroup("g1","g2","g3","g4","g5","g6","g7");;
# quotient out by relators to obtain fundamental group G of genus-3 surface
G:=F/[F.2*F.4*F.6*F.1*F.3*F.5*F.7,F.3*F.6*F.2*F.5*F.1*F.4*F.7];;
# set of generators of G and their inverses
gen:=[];;
gen[1]:=G.1;;gen[2]:=G.2;;gen[3]:=G.3;;gen[4]:=G.4;;gen[5]:=G.5;;gen[6]:=G.6;;gen[7]:=G.7;;
gen[8]:=G.1^-1;;gen[9]:=G.2^-1;;gen[10]:=G.3^-1;;gen[11]:=G.4^-1;;gen[12]:=G.5^-1;;gen[13]:=G.6^-1;;gen[14]:=G.7^-1;;
Print("done\n");

# define maximal index
n_max:=29;;
n_group:=3;
# compute all normal subgroups of index up to nmax
t_initial:=Runtime();;
Print("Computing first ",n_group," normal subgroups of index ",n_max,"... \n");
Graph_NSG:=LowIndexNormalSubgroupsSearchForIndex(G,n_max,n_group);;
List_NSG:=List(Graph_NSG);;
t_final:=Runtime();;
delta_time:=t_final-t_initial;
Print("Finished, it took ",StringTime(delta_time)," \n");

# number of normal subgroups found for each index
nsg:=[];;
for k in [1..n_max] do
	nsg[k]:=0;;
od;
for k in [1..Length(List_NSG)] do
	nsg[Index(List_NSG[k])]:=nsg[Index(List_NSG[k])]+1;
od;
for k in [1..n_max] do
	Print("Index ",k,": ",nsg[k]," normal subgroups\n");
od;
workspace_file:=Concatenation("./Workspaces/Index",String(n_max),".bin");;
Print("Saving workspace on ",workspace_file,"\n");
SaveWorkspace(workspace_file);;
