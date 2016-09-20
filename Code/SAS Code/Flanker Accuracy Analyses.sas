%let projectpath = C:\Documents and Settings\Michael Hallquist\My Documents\Online_Documents\Lab\Dissertation\Analyses;
%let datapath = &projectPath\Data;
*%let datapath = E:\Online_Documents\Lab\Dissertation\Analyses\Data;
libname diss "&datapath";


*determine frequency of incorrect responses and whether these loaded on incongruent and block type;
proc freq data = diss.flanker;
	tables CongruentBlock * Incongruent *  TrialSlide_ACC /chisq plcorr cellchi2;
run;

*accuracy does depend on trial type;
proc freq data = diss.flanker;
	tables Incongruent * TrialSlide_ACC /chisq plcorr cellchi2;
run;

*accuracy does depend on block type (more errors in incongruent blocks... but qualified by interaction?;
proc freq data = diss.flanker;
	tables CongruentBlock * TrialSlide_ACC /chisq plcorr cellchi2;
run;


proc logistic data = diss.flanker;
	model TrialSlide_ACC (event='0') = CongruentBlock | Incongruent;
quit;

*modeling event 0 means that we are modeling errors/inaccuracies;
data temp;
	set diss.flanker;
	int = CongruentBlock * Incongruent;
	*if subject in (36, 55, 83) then delete;
	if subject in (7, 14, 36, 38, 47, 55, 62, 83) then delete;
run;

proc logistic data = temp;
	*class CongruentBlock Incongruent;
	model TrialSlide_ACC (event='0') = CongruentBlock Incongruent int /expest;
	*exact 't1' Incongruent /estimate = odds;
quit;

proc glimmix data = temp;
	class CongruentBlock Incongruent;
	model TrialSlide_ACC = CongruentBlock Incongruent /solution dist=binary link=logit oddsratio;
	random intercept /sub=subject;
	nloptions technique=nrridg;
	*lsmeans Incongruent CongruentBlock /oddsratio;
run;

proc glimmix data = temp;
	class CongruentBlock Incongruent;
	model TrialSlide_ACC = CongruentBlock Incongruent /solution dist=binary link=logit oddsratio;
	random _residual_ /sub=subject type=ar(1);
	nloptions technique=nrridg;
	*lsmeans Incongruent CongruentBlock /oddsratio;
run;

proc glimmix data = temp;
	class CongruentBlock Incongruent;
	model TrialSlide_ACC = CongruentBlock Incongruent /solution dist=binary link=logit oddsratio;
	random _residual_ /sub=subject type=vc;
	nloptions technique=nrridg;
	*lsmeans Incongruent CongruentBlock /oddsratio;
run;

*compound symmetry makes the most sense;
*FINAL GLMM BELOW;
ods html;
ods graphics on;
proc glimmix data = temp plots=studentpanel(unpack);
	class CongruentBlock Incongruent;
	model TrialSlide_ACC = CongruentBlock | Incongruent /solution dist=binary link=logit oddsratio;
	random _residual_ /sub=subject type=cs;
	nloptions technique=nrridg;
	lsmeans Incongruent*CongruentBlock /odds oddsratio diff plots=meanplot(sliceby=congruentblock join);
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
	lsmestimate Incongruent -1 1 /exp cl;  *these estimate ORs > 1, easier to interp;
	lsmestimate CongruentBlock -1 1 /exp cl;

run;
ods graphics off;
ods html close;

proc sort data = glimOut;
	by subject;
run;

data residCheck;
	set glimOut;
	where student > 3;
run;

proc freq data = residCheck;
	tables subject;
run;

proc print data = glimOut;
	var subject student resid predicted trialslide_acc;
	where student > 4;
run;

proc print data = glimOut (firstobs=18000);
	var subject student resid predicted trialslide_acc;
run;



proc logistic data = diss.flankerScores;
	*class CongruentBlock Incongruent;
	model TrialSlide_ACC (event='0') = CongruentBlock Incongruent /expest;
quit;

proc freq data = diss.flankerScores (where=(TrialSlide_ACC=0));
	tables Incongruent /chisq plcorr cellchi2;
run;


proc sort data = diss.flankerScores;
	by subject;
run;


%inc "&projectPath\CorrWithMacro.sas";
*correlations of accuracy scores;
%corrWith(diss.allcombined, C_Inc_PctInaccurate I_Inc_PctInaccurate I_PctInaccurate PctInaccurate, SNAPFactors);
%corrWith(diss.allcombined, C_Inc_PctInaccurate I_Inc_PctInaccurate I_PctInaccurate C_PctInaccurate PctInaccurate, SNAPFacets);
%corrWith(diss.allcombined, C_Inc_PctInaccurate I_Inc_PctInaccurate I_PctInaccurate PctInaccurate, MPQFactors);
%corrWith(diss.allcombined, C_Inc_PctInaccurate I_Inc_PctInaccurate I_PctInaccurate PctInaccurate, MPQFacets);
%corrWith(diss.allcombined, C_Inc_PctInaccurate I_Inc_PctInaccurate I_PctInaccurate PctInaccurate, ATQ);


%corrWith(diss.allcombined, C_Inc_PctInaccurate I_Inc_PctInaccurate Inc_PctInaccurate PctInaccurate, SNAPFactors);
%corrWith(diss.allcombined, I_Con_PctInaccurate I_Inc_PctInaccurate I_PctInaccurate C_PctInaccurate PctInaccurate, SNAPFactors);

%corrWith(diss.allcombined, ATQ_ECTotal, MPQFactors);

%corrWith(diss.allcombined, Inc_PctInaccurate, IMPUL PROPER HDWK);


*start afresh;
ods html;
*%corrWith(diss.allcombined, C_Inc_PctInaccurate I_Inc_PctInaccurate Inc_PctInaccurate, SNAPFacets, rank=true, partial=Con_PctInaccurate);
%corrWith(diss.allcombined, C_I_Error I_I_Error I_Error, SNAPFacets, partial=C_Error);
ods html close;

ods html;
%corrWith(diss.allcombined, C_Inc_PctInaccurate I_Inc_PctInaccurate Inc_PctInaccurate, SNAPFacets, rank=true);
%corrWith(diss.allcombined, C_I_Error I_I_Error I_Error, SNAPFacets, rank=true);
ods html close;
%corrWith(diss.allcombined, C_Inc_PctInaccurate I_Inc_PctInaccurate Inc_PctInaccurate, SNAPFacets, rank=true);

%corrWith(diss.allcombined, accDiff, SNAPFacets);
%corrWith(diss.allcombined, Inc_PctInaccurate, SNAPFacets, partial=Con_PctInaccurate);

%corrWith(diss.allcombined, Inc_PctInaccurate, SNAPFacets, partial=Con_PctInaccurate);


proc corr data = diss.allcombined;
	var disinh disinhp;
run;

proc corr data = diss.allcombined;
	title 'intercorrelations among accuracy measures';
	var C_Con_PctInaccurate C_Inc_PctInaccurate I_Con_PctInaccurate I_Inc_PctInaccurate
		Inc_PctInaccurate Con_PctInaccurate C_PctInaccurate I_PctInaccurate PctInaccurate;
run;

proc corr data = diss.allcombined;
	var Inc_PctInaccurate disinh;
	partial Con_PctInaccurate;
run;

proc corr data = diss.allcombined;
	var Inc_PctInaccurate Con_PctInaccurate;
run;


proc corr data = diss.allcombined;
	var I_PctInaccurate disinh;
	partial C_PctInaccurate;
run;
