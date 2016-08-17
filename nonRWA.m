(* ::Package:: *)

BeginPackage["levelcrossing`"]

element::usage = "element=element[nmax,n,q]
	Takes n and q (and nmax), and gives proper element in my convention.
"
\[Psi]::usage = "barevector=\[Psi][nmax,tmax,n,q]
	Takes n and q (and nmax and tmax), and gives proper bare state in my convention.
"

Diagonalize::usage = "{Eigenvectors,Eigenenergies} = Diagonalize[\[Omega]q,\[Eta],\[Omega]r,g,g_num,ng,tmax,nmax]
	Takes system parameters and gives properly sorted eigenstates and eigenvectos.
	The qubit is a transmon, and is capacitively (charge-charge) coupled to the resonator.
	The coupling between qubit and resoantor is RWA, or in other words is excitation preserving.

		\[Omega]q: Qubit frequency.
		\[Eta]: Qubit anharmonicity.
		\[Omega]r: Resonator frequency.
		g: Coupling strenght between resonator and qubit, defined between |0,1> and |1,0> (|qubit-resonator>).
		g_num: Numerical g. This can be used to reduce the nmax in the system for faster calculations. This way any 'simulation n' in the code would correspond to physical photon number of n(g_num/g)^2
		ng: Background charge of transmon. WARNING: This value should not be of the form k/2. Always use non-excat values such as 0.4999.
		tmax: Maximum number of levels for the transmon. For typical parameters 9-10 is enough. WARNING: This value must be an integer.
		nmax: Maximum number of levels for the resonator. WARNING: This value must be an integer.

		Eigenvectors: A tmax*nmax by tmax*nmax matrix, including the eigenvectors. Each eigenvector is a column in this matrix, hence U|barebasis_element> = |eigenbasis_element>.
		Eigenenergies: A tmax*nmax list (column vector), including the eigenenergies.

	NAVIGATING THE RESULT
		Refer to GetEigen[] and GetBare[] functions manuals to find your desired eigenstate/eigenenergy and compare it to the bare state/energy
"

GetBare::usage = "{Barestate,Bareenergy}=GetBare[\[Omega]q,\[Eta],\[Omega]r,ng,tmax,nmax,q,n]
	Takes System Parameters and gives the the bare state/energy that belongs to qubit state ''0 <= q <= tmax-1'' and resonator state ''0 <= n <= nmax-1''

		\[Omega]q: Qubit frequency.
		\[Eta]: Qubit anharmonicity.
		\[Omega]r: Resonator frequency.
		ng: Background charge of transmon. WARNING: This value should not be of the form k/2. Always use non-excat values such as 0.4999.
		tmax: Maximum number of levels for the transmon. WARNING: This value must be an integer.
		nmax: Maximum number of levels for the resonator. WARNING: This value must be an integer.
		q: The qubit state that you want the eigenvalues for. This must be an integer in the range ''0 <= q <= tmax-1''.
		n: The resonator state that you want the eigenvalues for. This must be an integer in the range ''0 <= n <= nmax-1''.

		Barestate: The bare state belonging to qubit state q and resonator state n.
		Bareenergy: The bare energy belonging to qubit state q and resonator state n.
"

GetEigen::usage = "{Eigenstate,Eigenenergy}=GetEigen[U,E,tmax,nmax,q,n]
	Takes the ORDERED eigenstate and eigenenergies, and gives the the eigenenergy/eigenvector that belongs to qubit state ''0 <= q <= tmax-1'' and resonator state ''0 <= n <= nmax-1''

		U: Ordered matrix including all eigenvectors. Essentially output of Diagonalize.
		E: Ordered list including all eigenvectors. Essentially output of Diagonalize.
		tmax: Maximum number of levels for the transmon. WARNING: This value must be an integer.
		nmax: Maximum number of levels for the resonator. WARNING: This value must be an integer.
		q: The qubit state that you want the eigenvalues for. This must be an integer in the range ''0 <= q <= tmax-1''.
		n: The resonator state that you want the eigenvalues for. This must be an integer in the range ''0 <= n <= nmax-1''.

		Eigenstate: The eigenstate belonging to qubit state q and resonator state n.
		Eigenenergy: The eigenenergy belonging to qubit state q and resonator state n.
"
GetFan::usage = "{E0-vs-n,E1-vs-n,...Etmaxminus1-vs-n}=GetFan[E,tmax,nmax,\[Omega]r]
	Takes the ORDERED eigenenergies, and gives the the fan diagram of the system. This means the parameters should be consistent with the ones used in Diagonalize[]

		E: Ordered list including all eigenvectors. Essentially output of Diagonalize.
		tmax: Maximum number of levels for the transmon. WARNING: This value must be an integer.
		nmax: Maximum number of levels for the resonator. WARNING: This value must be an integer.
		\[Omega]r: Resonator frequency.

		Eq-vs-n: Energy in the fan diagram versus (simulation) photon number n for qubit q. To convert to physical n, multiply simulation n by (g_num/g)^2
"

GetCross::usage = "{n_0,n_1,...,n_tmaxminus1}=GetCross[E,tmax,nmax,\[Omega]r,initial_state,resonator_shift]
	Takes the ORDERED eigenenergies, and gives the the fan diagram of the system.

		E: Ordered list including all eigenvectors. Essentially output of Diagonalize.
		tmax: Maximum number of levels for the transmon. WARNING: This value must be an integer.
		nmax: Maximum number of levels for the resonator. WARNING: This value must be an integer.
		\[Omega]r: Resonator frequency.
		initial_state: The initial state that you want to find the level crossings from. WARNING: This must be an integer in the range ''0 <= q <= tmax-1''.
		resonator_shift: The number of resonator energies to shift for the level crossing. WARNING: This must be an integer, either 1 or 2.

		n_q: (simulation) photon number n at which there is a crossing between initial_state and the qubit state q. n_q=0 means there is no level crossing. Obviously, always n_q=0 for q<=initial_state.
			To convert to physical n at crossing, multiply n_q by (gnum/g)^2.
"

GetEffg::usage = "{{coh_geff_0,coh_geff_1,...,coh_geff_tmaxminus1},{incoh_geff_0,incoh_geff_1,...,incoh_geff_tmaxminus1}}=GetCross[Veigen,crossing_list,\[Omega]q,\[Eta],\[Omega]r,g,gnum,ng,tmax,nmax,initial_state,resonator_shift,broken_symmetry]
	Takes the eigenstates and list of n at which crossing happens, gives the effective coupling at the crossing point

		Veigen: The matrix of ordered eigenstates. Essentially the output of Diagonalize[].
		crossing_list: List of n at which crossing happens. This is the output of GetCross[].
		\[Omega]q: Qubit frequency.
		\[Eta]: Qubit anharmonicity.
		\[Omega]r: Resonator frequency.
		g: Coupling strenght between resonator and qubit, defined between |0,1> and |1,0> (|qubit-resonator>).
		g_num: Numerical g. This can be used to reduce the nmax in the system for faster calculations. This way any 'simulation n' in the code would correspond to physical photon number of n(g_num/g)^2
		tmax: Maximum number of levels for the transmon. WARNING: This value must be an integer.
		nmax: Maximum number of levels for the resonator. WARNING: This value must be an integer.
		\[Omega]r: Resonator frequency.
		initial_state: The initial state that the crossing starts from. This must be an integer in the range ''0 <= q <= tmax-1''.
		resonator_shift: The number of resonator energies to shift for the level crossing. This must be an integer, either 1 or 2.
		broken_symmetry: The breakage of selection rule. The program sets: \!\(\*SubscriptBox[\(g\), \(k, k + 2\)]\)=broken_symmetry*\!\(\*SqrtBox[\(\((k + 1)\) \((k + 2)\)\)]\).

		coh_geff_q: Coherent effective coupling normalized by g_num, for corssing between initial_state and qubit state q. Value of 0 means no crossing for that state.
		incoh_geff_q: Incoherent effective coupling normalized by g_num, for corssing between initial_state and qubit state q. Value of 0 means no crossing for that state.
"

Begin["Private`"]

element[nmax_,n_,q_]:=q nmax+(n+1);

\[Psi][nmax_,tmax_,n_, q_]:=SparseArray[{i_} /; i == (q nmax + (n + 1)) -> 1, {nmax*tmax}];(* State for resonator state n and qubit state q *)

GetBare[in\[Omega]q_,in\[Eta]_,in\[Omega]r_,inng_,intmax_,innmax_,inq_,inn_]:=
Block[{vector,energy,f,Et,r,Ec,\[Omega]t,\[Omega]q=in\[Omega]q,\[Eta]=in\[Eta],\[Omega]r=in\[Omega]r,ng=inng,tmax=intmax,nmax=innmax},
	f[k_]:=k+1-Mod[k+1,2]+2ng (-1)^(k-((Sign[ng]-1)/2));(*function that sorts eigenenergies of transmon (cosine potential) *)
	Et[k_,rr_,Ecc_]:=Ecc MathieuCharacteristicA[f[k],-1/2 rr];(*transmon 0\[LessEqual]kth eigenenergy. Ec is the capacitance energy, Ej is the juction energy, and r=EJ/Ec*)
	{r,Ec}=Block[{foundr,foundEc,result},
	result=FindRoot[{(Et[1,foundr,foundEc]-Et[0,foundr,foundEc])==\[Omega]q,((Et[1,foundr,foundEc]-Et[0,foundr,foundEc])-((Et[2,foundr,foundEc]-Et[1,foundr,foundEc])))==\[Eta]},{{foundr,100},{foundEc,0.2}}];
	{foundr,foundEc}/.result];(*this block find the correspinding r and Ec, given values of \[Omega]q and \[Eta] anf ng*)
	\[Omega]t=Table[Et[i,r,Ec]-Et[0,r,Ec],{i,0,tmax,1}];(*list of transmon energies*)

	vector=SparseArray[{i_}/;i==(inq innmax+(inn+1))->1,{innmax*intmax}];(* State for resonator state n and qubit state q *)
	energy=\[Omega]t[[inq+1]]+inn \[Omega]r;
{vector,energy}];

GetEigen[U_,E_,intmax_,innmax_,inq_,inn_]:=
Block[{nmax=innmax,vector,energy},
	energy=E[[element[nmax,inn,inq]]];
	vector=U\[Transpose][[element[nmax,inn,inq]]];
{vector,energy}];

Diagonalize[in\[Omega]q_,in\[Eta]_,in\[Omega]r_,ing_,ingnum_,inng_,intmax_,innmax_]:=
Block[{\[Omega]q=in\[Omega]q,\[Eta]=in\[Eta],\[Omega]r=in\[Omega]r,g=ing,gnum=ingnum,ng=inng,tmax=intmax,nmax=innmax,r,Ec,f,\[Psi]t,Et,\[Omega]t,gnorm,gkkp1,elem,Hr,Hq,HRWA,U,Eeigen},

f[k_]:=k+1-Mod[k+1,2]+2ng (-1)^(k-((Sign[ng]-1)/2));(*function that sorts eigenenergies of transmon (cosine potential) *)
Et[k_,rr_,Ecc_]:=Ecc MathieuCharacteristicA[f[k],-1/2 rr];(*transmon 0\[LessEqual]kth eigenenergy. Ec is the capacitance energy, Ej is the juction energy, and r=EJ/Ec*)
{r,Ec}=Block[{foundr,foundEc,result},
result=FindRoot[{(Et[1,foundr,foundEc]-Et[0,foundr,foundEc])==\[Omega]q,((Et[1,foundr,foundEc]-Et[0,foundr,foundEc])-((Et[2,foundr,foundEc]-Et[1,foundr,foundEc])))==\[Eta]},{{foundr,100},{foundEc,0.2}}];
{foundr,foundEc}/.result];(*this block find the correspinding r and Ec, given values of \[Omega]q and \[Eta] anf ng*)
\[Psi]t[k_,\[CurlyPhi]_]:=E^(I ng \[CurlyPhi])/Sqrt[2 \[Pi]] (MathieuC[MathieuCharacteristicA[f[k],-(r/2)],-(r/2),\[CurlyPhi]/2]+I (-1)^(k+1) MathieuS[MathieuCharacteristicA[f[k],-(r/2)],-(r/2),\[CurlyPhi]/2]);(*transmon 0\[LessEqual]kth eigenstate as a function of phase \[CurlyPhi]*)

\[Omega]t=Table[Et[i,r,Ec]-Et[0,r,Ec],{i,0,tmax,1}];(*list of transmon energies*)
gnorm=Chop[NIntegrate[\[Psi]t[0,\[CurlyPhi]]\[Conjugate] 1/I D[\[Psi]t[1,\[CurlyPhi]],\[CurlyPhi]],{\[CurlyPhi],-\[Pi],\[Pi]}]];(*normalization factor for coupling*)
gkkp1=Table[Abs[gnum/gnorm Chop[NIntegrate[\[Psi]t[i,\[CurlyPhi]]\[Conjugate] 1/I D[\[Psi]t[i+1,\[CurlyPhi]],\[CurlyPhi]],{\[CurlyPhi],-\[Pi],\[Pi]}]]],{i,0,tmax-2,1}];(*list of Subscript[g, k,k+1] couplings*)

elem[n_,q_]:=q nmax+(n+1);(*gives proper vector element for resonator state n\[Element][0,nmax-1] and qubit state q\[Element][0,tmax-1]. This sets my convention for naming my states*)

Hr=\[Omega]r SparseArray[Flatten[Table[Table[{elem[n,q],elem[n,q]}->n,{n,0,nmax-1}],{q,0,tmax-1}]]];(*Resonator Hamiltonian*)
Hq=SparseArray[Flatten[Table[Table[{elem[n,q],elem[n,q]}->\[Omega]t[[q+1]],{n,0,nmax-1}],{q,0,tmax-1}]]];(*Qubit Hamiltonian*)
HRWA=SparseArray[(Flatten[Table[Table[{elem[n,q],elem[n-1,q+1]}->gkkp1[[q+1]]Sqrt[n],{n,1,nmax-1}],{q,0,tmax-2}]]),{nmax*tmax,nmax*tmax}]+(SparseArray[(Flatten[Table[Table[{elem[n,q],elem[n-1,q+1]}->gkkp1[[q+1]]Sqrt[n],{n,1,nmax-1}],{q,0,tmax-2}]]),{nmax*tmax,nmax*tmax}])\[ConjugateTranspose];(*RWA interaction Hamiltonian*)

(*The block below diagonalizes the system*)
Block[{Uul={},eigenenergy={},bareenergy={},Ud={},eigenenergyordered={},intlist={},eigenpositions={},barepositions={},eigenorder={},bareorder={},ref,system,amps,sgn,ladder},
	system=Quiet[Eigensystem[N[Hr+Hq+HRWA]]];
	Uul=SparseArray[Chop[system[[2]],10^-7]];(*Unordered unitary diagonalizer matrix, Uul\[ConjugateTranspose].H.Uul=Subscript[H, D]*)
	eigenenergy=system[[1]];(*Unordered eigenenergy*)
	bareenergy=Normal[Diagonal[N[Hr+Hq]]];
	Ud=SparseArray[{{1,1}->1,{nmax*tmax,nmax*tmax}->1}];
	eigenenergyordered=SparseArray[{{1}->bareenergy[[1]],{nmax*tmax}->bareenergy[[nmax*tmax]]}];
	ladder[i_]:=Piecewise[{{i,i<= tmax-1},{tmax-1,tmax<=i<= nmax-1}}];(*function to take care of variable interacting subspace*)
	intlist=Table[Table[elem[n-q,q],{q,0,ladder[n]}],{n,1,nmax-1}]~Join~Table[Table[elem[(nmax-1)-(q-qi),q],{q,qi,tmax-1}],{qi,1,tmax-2}];(*list of interacting states*)
	eigenpositions=Table[Union[Flatten[(Uul[[;;,#]]["NonzeroPositions"]&/@intlist[[i]])]],{i,1,Length[intlist]}];(*position of the eigenstates belonging to the same interacting subspace*)
	barepositions=intlist;(*position of the bare states belonging to the same inbteracting subspace*)
	eigenorder=Ordering[(eigenenergy[[eigenpositions[[#]]]])]&/@Range[Length[eigenpositions]];(*order of eigenstate energies*)
	bareorder=Ordering[bareenergy[[barepositions[[#]]]]]&/@Range[Length[eigenpositions]];(*order of bare state energies*)
	(Ud[[barepositions[[#]]]]=Uul[[(Permute[eigenpositions[[#]],FindPermutation[eigenorder[[#]],bareorder[[#]]]])]])&/@Range[Length[eigenpositions]];(*reordering eigenstates*)
	(eigenenergyordered[[barepositions[[#]]]]=eigenenergy[[(Permute[eigenpositions[[#]],FindPermutation[eigenorder[[#]],bareorder[[#]]]])]])&/@Range[Length[eigenpositions]];(*reordering eigenenergies*)
	Do[(*correcting the eigenstate amplitude signs*)
		ref=Sign[Ud[[elem[0,q],elem[0,q]]]];(*reference sign, chosen to be the same as the first element*)
		(Ud[[elem[#,q]]]=ref*Sign[Ud[[elem[#,q],elem[#,q]]]]*Ud[[elem[#,q]]])&/@Range[nmax-1];(*correcting the signs for all n*)
	,{q,0,tmax-1}](*doing for all qubit states*);
	Do[(*correcting the sign when fist element changes sign, by detecting when this sign change happens*)
		If[q!=tmax-1,amps=Ud[[elem[#,q],elem[#-1,q+1]]]&/@Range[1,nmax-1],amps=Ud[[elem[#,q],elem[#+2,q-2]]]&/@Range[0,nmax-2]];
		ref=Sign[Ud[[elem[0,q],elem[0,q]]]];
		sgn=ref*Sign[amps];
		(Ud[[elem[#,q]]]=sgn[[#]]*Ud[[elem[#,q]]])&/@Range[nmax-1];
	,{q,0,tmax-1}](*doing for all qubit states*);
	{SparseArray[Ud\[ConjugateTranspose]],eigenenergyordered}]
];

GetFan[E_,intmax_,innmax_,in\[Omega]r_]:=
Block[{Eeigen=E,tmax=intmax,nmax=innmax,\[Omega]r=in\[Omega]r,fan},
	fan=Table[{n,Eeigen[[element[nmax,n-q,q]]]-n*\[Omega]r},{q,0,tmax-1},{n,q,nmax-5}];
fan];

GetCross[E_,intmax_,innmax_,in\[Omega]r_,ininit_,inshift_]:=
Block[{Eeigen=E,tmax=intmax,nmax=innmax,\[Omega]r=in\[Omega]r,init=ininit,shift=inshift,fan,difflist,ncross,poscross},
	fan = GetFan[Eeigen, tmax, nmax, \[Omega]r];
	(*|init,n+q> --> |init+q,n+shift> + shift*\[Omega]r*)
	difflist = Table[
					Table[
						{fan[[init + 1, (n + q) + 1, 1]], fan[[init + 1, (n + q) + 1, 2]]-fan[[(init + q) + 1, (n + shift) + 1, 2]] - shift \[Omega]r}
					,{n, 0, Length[fan[[init + q + 1]]] - shift - 1}]
				,{q,  1, tmax - init - 1}];
	difflist = Abs[difflist];(*find a list of substraction of fan curves*)
	ncross = Table[
				poscross = Position[difflist[[p]]\[Transpose][[2]],Min[difflist[[p]]\[Transpose][[2]]]][[1, 1]];(*the crossing happens when the difference between curves becomes zero*)
				If[poscross == Length[difflist[[p]]] || poscross == 1, 0 , difflist[[p]]\[Transpose][[1,poscross]]](*make sure we are throwing out the detection which are at the very end or very begining*)
			,{p, 1, Length[difflist]}];
Join[ConstantArray[0, init + 1], ncross]](*Prepend the nq=0 for q\[LessEqual] initial_state and give as output*)

GetEffg[U_,inncross_,in\[Omega]q_,in\[Eta]_,in\[Omega]r_,ing_,ingnum_,inng_,intmax_,innmax_,ininit_,inshift_,inviolation02_]:=
Block[{ncross=inncross,\[Omega]q=in\[Omega]q,\[Eta]=in\[Eta],\[Omega]r=in\[Omega]r,g=ing,gnum=ingnum,ng=inng,tmax=intmax,nmax=innmax,init=ininit,shift=inshift,violation02=inviolation02,f,Et,r,Ec,\[Psi]t,\[Omega]t,gnorm,gkkp1,gkkp3,gkkp2,data,nx,H1list,H3list,H2list},
	f[k_]:=k+1-Mod[k+1,2]+2ng (-1)^(k-((Sign[ng]-1)/2));(*function that sorts eigenenergies of transmon (cosine potential)*)
	Et[k_,rr_,Ecc_]:=Ecc MathieuCharacteristicA[f[k],-1/2 rr];(*transmon 0\[LessEqual]kth eigenenergy.Ec is the capacitance energy,Ej is the juction energy,and r=EJ/Ec*)
	{r,Ec}=Block[{foundr,foundEc,result},result=FindRoot[{(Et[1,foundr,foundEc]-Et[0,foundr,foundEc])==\[Omega]q,((Et[1,foundr,foundEc]-Et[0,foundr,foundEc])-((Et[2,foundr,foundEc]-Et[1,foundr,foundEc])))==\[Eta]},{{foundr,100},{foundEc,0.2}}];
		{foundr,foundEc}/.result];(*this block find the correspinding r and Ec,given values of \[Omega]q and \[Eta] anf ng*)

	\[Psi]t[k_,\[CurlyPhi]_]:=E^(I ng \[CurlyPhi])/Sqrt[2 \[Pi]] (MathieuC[MathieuCharacteristicA[f[k],-(r/2)],-(r/2),\[CurlyPhi]/2]+I (-1)^(k+1) MathieuS[MathieuCharacteristicA[f[k],-(r/2)],-(r/2),\[CurlyPhi]/2]);(*transmon 0\[LessEqual]kth eigenstate as a function of phase \[CurlyPhi]*)
	\[Omega]t=Table[Et[i,r,Ec]-Et[0,r,Ec],{i,0,tmax,1}];(*list of transmon energies*)
	gnorm=Chop[NIntegrate[\[Psi]t[0,\[CurlyPhi]]\[Conjugate] 1/I D[\[Psi]t[1,\[CurlyPhi]],\[CurlyPhi]],{\[CurlyPhi],-\[Pi],\[Pi]}]];(*normalization factor for coupling*)
	gkkp1=Table[Abs[gnum/gnorm Chop[NIntegrate[\[Psi]t[i,\[CurlyPhi]]\[Conjugate] 1/I D[\[Psi]t[i+1,\[CurlyPhi]],\[CurlyPhi]],{\[CurlyPhi],-\[Pi],\[Pi]}]]],{i,0,tmax-2,1}];(*list of Subscript[g, k,k+1] couplings*)
	gkkp3=Table[Abs[ (gnum/gnorm)Chop[NIntegrate[\[Psi]t[i,\[CurlyPhi]]\[Conjugate] 1/I D[\[Psi]t[i+3,\[CurlyPhi]],\[CurlyPhi]],{\[CurlyPhi],-\[Pi],\[Pi]}]]],{i,0,(tmax-4)}];(*list of Subscript[g, k,k+3] couplings*)
	gkkp2=Table[Abs[ gnum violation02 Sqrt[(i+1)(i+2)]],{i,0,(tmax-3)}];(*list of Subscript[g, k,k+2] couplings*)

Switch[shift,1,(*When shift is for 1 RWA strip*)
	data=Reap[
			Do[
				If[ncross[[q+1]]==0,Sow[{0,0}],
					nx=ncross[[q+1]];
					H1list=Table[(\[Psi][nmax,tmax,init+nx-l,l].(U.\[Psi][nmax,tmax,nx,init]))gkkp1[[l+1]]Sqrt[init+nx-l+1](\[Psi][nmax,tmax,init+nx+shift-l-1,l+1].(U.\[Psi][nmax,tmax,init+nx+shift-q,q]))\[Conjugate],{l,0,tmax-2}];
					H3list=Table[(\[Psi][nmax,tmax,init+nx-l,l].(U.\[Psi][nmax,tmax,nx,init]))gkkp3[[l+1]]Sqrt[init+nx-l](\[Psi][nmax,tmax,init+nx+shift-l-3,l+3].(U.\[Psi][nmax,tmax,init+nx+shift-q,q]))\[Conjugate],{l,0,tmax-4}];
					H2list=Table[(\[Psi][nmax,tmax,init+nx-l,l].(U.\[Psi][nmax,tmax,nx,init]))gkkp2[[l+1]]Sqrt[init+nx-l](\[Psi][nmax,tmax,init+nx+shift-l-2,l+2].(U.\[Psi][nmax,tmax,init+nx+shift-q,q]))\[Conjugate],{l,0,tmax-3}];
					Sow[Abs[{Total@H2list,Norm[H2list]}/gnum]];
					]
			,{q,0,Length[ncross]-1}];
			]//Last//First;
,2,(*When shift is for 2 RWA strip*)
	data=Reap[
			Do[
				If[ncross[[q+1]]==0,Sow[{0,0}],
					nx=ncross[[q+1]];
					H1list=Table[(\[Psi][nmax,tmax,init+nx-l,l].(U.\[Psi][nmax,tmax,nx,init]))gkkp1[[l+1]]Sqrt[init+nx-l+1](\[Psi][nmax,tmax,init+nx+shift-l-1,l+1].(U.\[Psi][nmax,tmax,init+nx+shift-q,q]))\[Conjugate],{l,0,tmax-2}];
					H3list=Table[(\[Psi][nmax,tmax,init+nx-l,l].(U.\[Psi][nmax,tmax,nx,init]))gkkp3[[l+1]]Sqrt[init+nx-l](\[Psi][nmax,tmax,init+nx+shift-l-3,l+3].(U.\[Psi][nmax,tmax,init+nx+shift-q,q]))\[Conjugate],{l,0,tmax-4}];
					H2list=Table[(\[Psi][nmax,tmax,init+nx-l,l].(U.\[Psi][nmax,tmax,nx,init]))gkkp2[[l+1]]Sqrt[init+nx-l](\[Psi][nmax,tmax,init+nx+shift-l-2,l+2].(U.\[Psi][nmax,tmax,init+nx+shift-q,q]))\[Conjugate],{l,0,tmax-3}];
					Sow[Abs[{Total@H1list+Total@H3list,(Norm[H1list]^2+Norm[H3list]^2)^(1/2)}/gnum]];
					]
			,{q,0,Length[ncross]-1}];
			]//Last//First;
];
data\[Transpose]];

End[]
EndPackage[]
