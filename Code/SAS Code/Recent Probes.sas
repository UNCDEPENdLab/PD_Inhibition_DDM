%let projectpath = C:\Documents and Settings\Michael Hallquist\My Documents\Online_Documents\Lab\Dissertation\Analyses;
*%let projectpath = E:\Online_Documents\Lab\Dissertation\Analyses;
%let datapath = &projectPath\Data;
libname diss "&datapath";


proc means data = diss.RecentProbes;
	class condition;
	var ProbeDisplay_RT;
run;

proc univariate data = diss.RecentProbesScores plots;
	var RP_ACC;
	id subject;
run;

proc univariate data = diss.RecentProbesScores plots;
	var RP_RT;
	id subject;
run;

*proc freq data = diss.RecentProbes; *(where=(condition in ("Yfam0" "Yfam1" "Yfam1RI" "Yfam2")));
*	tables ProbeDisplay_ACC * condition /chisq CELLCHI2;
*run;

proc print data = diss.RecentProbesScores;
run;

proc contents data = diss.RecentProbesScores;
run;

proc glm data = diss.RecentProbesScores;
	model RP_Nfam0_RT RP_Nfam1_RT RP_Nfam1RI_RT RP_Nfam2_RT RP_Y_RT = /nouni;
	contrast 'nfam0 v nfam1' condition 1 -1 0 0 0;
	contrast 'n resp v y resp' condition -0.25 -0.25 -0.25 -0.25 1;
	repeated condition 5 /summary; *contrast(1) /printm summary;
quit;


proc transpose data = diss.RecentProbesScores out=repMeasures;
	by subject;
	var RP_Nfam0_RT RP_Nfam1_RT RP_Nfam1RI_RT RP_Nfam2_RT RP_Y_RT;
run;

proc print data = repMeasures;
run;

proc print data = repMeasures (obs=10);
run;

*un is "unstructured", so the error covariances are uniquely estimated for each pair of levels;
proc mixed data = repMeasures covtest;
	class _NAME_ subject;
	model col1 = _NAME_;
	*random intercept /sub=subject;
	repeated /type = un subject=subject sscp rcorr;
	lsmeans _NAME_ /adjust=tukey pdiff;
	estimate 'unfamiliar versus familiar' _NAME_ -1 0.333 0.333 0.334 0 /E;
	estimate 'response interference v. familiar' _NAME_ 0 1 -0.5 -0.5 0 /E;
quit;

*vc is "variance components", basically the error covariances among levels are fixed at 0;
proc mixed data = repMeasures covtest;
	class _NAME_ subject;
	model col1 = _NAME_;
	*random intercept /sub=subject;
	repeated /type =vc subject=subject rcorr;
	lsmeans _NAME_ /adjust=tukey pdiff;
quit;

*cs is "compound symmetry", basically estimating a single covariance among all pairs;
proc mixed data = repMeasures covtest;
	class _NAME_ subject;
	model col1 = _NAME_;
	*random intercept /sub=subject;
	repeated /type =cs subject=subject rcorr;
	lsmeans _NAME_ /adjust=tukey pdiff;
	estimate 'unfamiliar versus familiar' _NAME_ -1 0 0.5 0.5 0 /E;
quit;



proc mixed data = repMeasures covtest;
	class _NAME_ subject;
	model col1 = _NAME_;
	random intercept /sub=subject;
	*repeated /type = un subject=subject rcorr;
	lsmeans _NAME_ /adjust=tukey pdiff;
quit;

proc mixed data = repMeasures covtest;
	class _NAME_ subject;
	model col1 = _NAME_;
	random intercept /sub=subject;
	repeated /type = un subject=subject rcorr;
	lsmeans _NAME_ /adjust=tukey pdiff;
quit;

proc sort data = diss.questionnaires;
	by subject;
run;

proc sort data = diss.RecentProbesScores;
	by subject;
run;

data merged;
	merge diss.questionnaires diss.RecentProbesScores;
	by subject;
run;

proc corr data = merged;
	var RP_Nfam0_RT RP_Nfam1_RT RP_Nfam1RI_RT RP_Nfam2_RT RP_Y_RT;
	with postemp negtemp disinh mps_per mps_ner mps_cor;
run;

proc corr data = merged (where=(mps_inv=0));
	var RP_Nfam0_RT RP_Nfam1_RT RP_Nfam1RI_RT RP_Nfam2_RT RP_Y_RT;
	with mps_wbr mps_spr mps_acr mps_scr mps_alr mps_agr mps_clr mps_har atq_ec;
run;


proc corr data = merged (where=(mps_inv=0));
	var RP_Nfam0_RT RP_Nfam1_RT RP_Nfam1RI_RT RP_Nfam2_RT RP_Y_RT;
	with ATQ_InhibControl ATQ_InhibControl ATQ_ActControl;
run;

proc corr data = merged (where=(mps_inv=0));
	var RP_Nfam0_ACC RP_Nfam1_ACC RP_Nfam1RI_ACC RP_Nfam2_ACC RP_Y_ACC;
	with ATQ_InhibControl ATQ_InhibControl ATQ_ActControl;
run;


proc print data = accOut;
run;

proc corr data = merged;
	var RP_Nfam0_ACC RP_Nfam1_ACC RP_Nfam1RI_ACC RP_Nfam2_ACC RP_Y_ACC;
	with postemp negtemp disinh mps_per mps_ner mps_cor;
run;

proc univariate data = diss.RecentProbesScores plots;
	var RP_Nfam0_ACC RP_Nfam1_ACC RP_Nfam1RI_ACC RP_Nfam2_ACC RP_Y_ACC;
run;


data temp;
	set diss.RecentProbes;
	if subject in (40, 65, 100) then delete; *near 50% in terms of accuracy;
	if condition = "Y" then delete;
run;

proc means data = diss.RecentProbesScores;
	title 'response accuracy rates across participants';
	var RP_Y_ACC
		RP_Nfam0_ACC
		RP_Nfam1_ACC
		RP_Nfam1RI_ACC
		RP_Nfam2_ACC
		RP_ACC;
run;


proc corr data = diss.allcombined;
	title 'response accuracy rates across participants';
	var RP_Nfam1_ACC
		RP_Nfam1RI_ACC
		RP_Nfam2_ACC;
run;

proc corr data = diss.RecentProbesScores;
	title 'response accuracy rates across participants';
	var RP_Nfam1_SUM
		RP_Nfam1RI_SUM
		RP_Nfam2_SUM
		RP_unfamiliar_SUM
		RP_Y_SUM
;
run;



proc corr data = diss.allcombined;
	title 'response accuracy rates across participants';
	var RP_Nfam1_RT
		RP_Nfam1RI_RT
		RP_Nfam2_RT
		RP_familiar_RT
		RP_Y_RT
	;
run;

*FINAL GLMM FOR RESPONSE ACCURACY;
ods html;
ods graphics on;
proc glimmix data = temp plots=studentpanel(unpack);
	class condition;
	model ProbeDisplay_ACC = condition /solution dist=binary link=logit oddsratio;
	random _residual_ /sub=subject type=cs;
	nloptions technique=nrridg;
	lsmeans condition /odds oddsratio adjust=tukey diff plots=meanplot(join);
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
	lsmestimate condition -1 0.5 0 0.5 0 /exp e;
	lsmestimate condition 'familiar v all unfamiliar' -3 1 1 1 /exp e divisor=3;

run;
ods graphics off;
ods html close;

*normal model to aid interpretation;
ods html;
ods graphics on;
proc glimmix data = temp plots=studentpanel(unpack);
	class condition;
	model ProbeDisplay_ACC = condition /solution dist=normal link=identity;
	random _residual_ /sub=subject type=cs;
	nloptions technique=nrridg;
	lsmeans condition /diff adjust=tukey plots=meanplot(join);
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
*	lsmestimate Incongruent -1 1 /exp cl;  *these estimate ORs > 1, easier to interp;
*	lsmestimate CongruentBlock -1 1 /exp cl;
	lsmestimate condition -1 0.5 0 0.5 0 /e;
	lsmestimate condition 'familiar v all unfamiliar' -3 1 1 1 /e divisor=3;

run;
ods graphics off;
ods html close;

proc sort data = glimOut;
	by subject;
run;

data residCheck;
	set glimOut;
	where abs(student) > 2;
run;

proc freq data = residCheck;
	tables subject;
run;

proc print data = diss.RecentProbesScores;
	var subject RP_Nfam0_ACC RP_Nfam1_ACC RP_Nfam1RI_ACC RP_Nfam2_ACC RP_Y_ACC RP_ACC;
	where subject in (14, 33, 38, 40, 65, 87, 100);
run;

proc transpose data = diss.RecentProbesScores out=RP_ACC;
	by subject;
	var RP_Nfam0_ACC RP_Nfam1_ACC RP_Nfam1RI_ACC RP_Nfam2_ACC RP_Y_ACC;
run;

data RP_ACC;
	set RP_ACC;
	if subject in (40, 65, 100) then delete; *near 50% in terms of accuracy;
	*if _NAME_ = "RP_Y_ACC" then delete;
run;

proc mixed data = RP_ACC covtest;
	class _NAME_ subject;
	model col1 = _NAME_ /outp=accOut residual influence;
	*random intercept /sub=subject;
	repeated /type = cs subject=subject rcorr;
	lsmeans _NAME_ /adjust=tukey pdiff;
	*estimate "unfamiliar versus familiar negative trials" _NAME_ -3 1 1 1 0;
quit;

proc print data = RP_ACC;
	var subject;
run;


%inc "&projectPath\CorrWithMacro.sas";
%corrWith(diss.allcombined (where=(subject not in(40, 65, 100))), RP_Nfam0_RT RP_Nfam1_RT RP_Nfam1RI_RT RP_Nfam2_RT RP_Y_RT, SNAPFacets);

%corrWith(diss.allcombined, RP_Nfam0_RT RP_Nfam1_RT RP_Nfam1RI_RT RP_Nfam2_RT RP_Y_RT, SNAPFacets, rank=true);
%corrWith(diss.allcombined, RP_familiar_RT, SNAPFacets, rank=true);
%corrWith(diss.allcombined, RP_familiar_RT, SNAPFacets, partial=RP_Y_RT, rank=true);
%corrWith(diss.allcombined, RP_familiar_SUM, SNAPFacets, partial=RP_Y_SUM, rank=true);
%corrWith(diss.allcombined, RP_familiar_SUM, SNAPFacets, rank=true);

%corrWith(diss.allcombined, RP_Nfam0_RT RP_Nfam1_RT RP_Nfam1RI_RT RP_Nfam2_RT RP_Y_RT, SNAPFactors);
%corrWith(diss.allcombined (where=(subject not in(40, 65, 100))), RP_Familiar_ACC RP_Nfam0_ACC RP_Nfam1_ACC RP_Nfam1RI_ACC RP_Nfam2_ACC RP_Y_ACC, SNAPFacets, rank=true);
%corrWith(diss.merged, RP_Familiar_ACC RP_Nfam0_ACC RP_Nfam1_ACC RP_Nfam1RI_ACC RP_Nfam2_ACC RP_Y_ACC, SNAPFactors);

data temp;
	set diss.RecentProbes;
	if subject in (40, 65, 100) then delete; *near 50% in terms of accuracy;
	if ProbeDisplay_ACC = 0 Then ProbeDisplay_RT = .;
	if condition = "Y" then delete;
run;

*GLIMMIX FOR REACTION TIME;
ods html;
ods graphics on;
proc glimmix data = temp plots=studentpanel(unpack) ic=pq;
	title 'recent probes RT as a function of trial type';
	class condition subject;
	model ProbeDisplay_RT = condition /solution dist=normal link=identity;
	random _residual_ /sub=subject type=cs;
	nloptions technique=nrridg;
	lsmeans condition /diff adjust=tukey plots=meanplot(sliceby=condition join);
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
	lsmestimate condition 'unfamiliar v familiar' -1 0.5 0 0.5 /e;
	lsmestimate condition 'unfamiliar v all familiar' -3 1 1 1 /e divisor=3;
run;
ods graphics off;
ods html close;

proc means data = temp;
	class condition;
	var ProbeDisplay_RT;
run;
