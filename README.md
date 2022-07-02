# -ARL-method
solve GAP by ARL method in GAMS
*************** sets ***************
Sets
I /i1*i3/
J /j1*j5/

**************  parameters  ********

Parameters

w(i,j)
t(i)
;

w(i,j) = uniform(0,1);
t(i)=uniform(1,1.7);
Display
w
t
*************  Variables  **********
Binary Variable
x(i,j)
;

Free Variable
Z_ARL
;
*************  model  **************
Equations
obj
cons1

;

Scalars
mio
Landa;

obj.. Z_ARL =e= sum({i,j},w(i,j)*x(i,j)) +sum({j},(mio/2)*power(sum({i},x(i,j)-1),2))+ sum({j},Landa*(sum({i},x(i,j))-1));

cons1(i)..  sum({j},w(i,j)*x(i,j)) =l= t(i);

*********

Model GAP_ALR
/
obj
cons1
/
;

Options
RESLIM=100
OPTCR=0
QCP=CPLEX
;
************** Setting

Set iter /1*100/

Parameters
MRE /0.01/
RE
Phi /0.0001/
Landa /0/
mio / 0/
Result(iter,*)
FC 'feasibility checker'
convergency
UB /10000/
LB /0.0001/
;
convergency = NO;
**************


LOOP(iter$(NOT(convergency)),

Solve  GAP_ALR us MIQCP max Z_ARL;

Result(iter,'L') = Landa;
Result(iter,'M') = mio;
LB=Z_ARL.l;
Result(iter,'LB = Z_ARL(L)') = Z_ARL.l;

** FC and LB
IF(sum({j},Landa*(sum({i},x.l(i,j))-1)),

FC=NO;

ELSE
FC=YES;
UB = sum({i,j},w(i,j)*x.l(i,j)) ;
)
;
**

Result(iter,'FC') = FC;
Result(iter,'LB') = LB;

RE = (UB - LB)/LB;

Result(iter,'RE') = RE;

*

IF( RE <= MRE ,
convergency = YES;
);


** update L

Landa = max { 0 , Landa  + Phi*(UB-LB)*(sum({j},Landa*(sum({i},x.l(i,j))-1))/abs(sum({j},Landa*(sum({i},x.l(i,j))-1))))}


)
*End of Loop
;

Display
Result
;



Execute_Unload 'ALR_IP401'
