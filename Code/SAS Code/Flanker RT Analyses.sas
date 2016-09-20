%let projectpath = C:\Documents and Settings\Michael Hallquist\My Documents\Online_Documents\Lab\Dissertation\Analyses;
%let datapath = &projectPath\Data;
*%let datapath = E:\Online_Documents\Lab\Dissertation\Analyses\Data;
libname diss "&datapath";

*model differences between incongruent and congruent trials using one-way random-effects ANOVA;
data flanker;
	set diss.flanker;
	if subject in (7, 14, 36, 38, 47, 55, 62, 83) then delete;
	retain i;
	if _n_ = 1 then i = 0;
	*for inaccurate responses, remove that RT from the running;
	if TrialSlide_ACC = 0 then do; 
		i = i+1;
		TrialSlide_RT = .;
	end;
	put "num inaccurate: " i;
	
run;

proc mixed data = flanker;
	class subject Incongruent;
	model TrialSlide_RT = Incongruent /solution;
	random intercept /sub=subject;
run;

ods trace on;
*two-way random effects ANOVA;
proc mixed data = flanker;
	class subject Incongruent CongruentBlock;
	model TrialSlide_RT = Incongruent | CongruentBlock /solution;
	random intercept /sub=subject;
	lsmeans Incongruent*CongruentBlock /adjust=tukey pdiff;
	ods output LSMeans = flankerInteraction;
run;
ods trace off;

ods trace on;
*two-way random effects ANOVA;
proc mixed data = flanker;
	class subject Incongruent CongruentBlock;
	model TrialSlide_RT = Incongruent | CongruentBlock /solution;
	repeated /sub=subject;
	lsmeans Incongruent*CongruentBlock /adjust=tukey pdiff;
	ods output LSMeans = flankerInteraction;
run;
ods trace off;

ods html;
ods graphics on;
proc glimmix data = flanker plots=studentpanel(unpack);
	class subject Incongruent CongruentBlock;
	model TrialSlide_RT = Incongruent | CongruentBlock /solution dist=normal link=identity;
	random _residual_ /sub=subject type=cs;
	nloptions technique=nrridg;
	lsmeans Incongruent /diff pdiff adjust=tukey plots=meanplot(join);
	lsmeans CongruentBlock /diff pdiff adjust=tukey plots=meanplot(join);
	lsmeans Incongruent*CongruentBlock /diff slice=incongruent slicediff=incongruent adjust=tukey plots=meanplot(sliceby=incongruent join);
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
	*lsmestimate BlockType -1 1 /exp cl;  *these estimate ORs > 1, easier to interp;

run;
ods graphics off;
ods html close;


proc print data = flankerInteraction;
run;

*interaction plot;
proc plot data = flankerInteraction;
	plot Estimate * Incongruent = CongruentBlock;
quit;

*create person-level mean values for congruent and incongruent trials;
*BEGIN HARMONIC MEAN COMPUTATION;
proc sort data = flanker;
	by subject CongruentBlock Incongruent;
run;

proc transpose data = flanker out=transFlanker; *(rename=(col1=I_congruentMean col2=C_congruentMean col3 = I_incongruentMean col4 = C_incongruentMean));
	var TrialSlide_RT;
	by subject CongruentBlock Incongruent;
run;

proc sort data = flankerTrim;
	by subject CongruentBlock Incongruent;
run;
proc transpose data = flankerTrim out=transTrim;
	var TrialSlide_RT;
	by subject CongruentBlock Incongruent;
run;

proc print data = transFlanker (where=(subject=105));
run;

data transFlanker;
	set transFlanker;
	harmean = harmean(of col1-col84);
	drop _name_ _label_;
run;


data transTrim;
	set transTrim;
	harmean = harmean(of col1-col84);
	drop _name_ _label_;
run;

*This mixed model compares harmonic means across block type and trial type;
*repeated measures approach, as opposed to random effects above;
proc mixed data = transFlanker covtest;
	class CongruentBlock Incongruent;
	model harmean = CongruentBlock | Incongruent /ddfm=kr solution;
	*random intercept /sub=subject;
	repeated /type = un subject=subject sscp rcorr;
	lsmeans CongruentBlock*Incongruent /adjust=tukey pdiff;
	ods output LSMeans = flankerInteraction;
	estimate 'con v incon' Incongruent 1 -1;
quit;

proc mixed data = transFlanker covtest;
	class Incongruent;
	model harmean = Incongruent /ddfm=kr solution;
	*random intercept /sub=subject;
	repeated /type = un subject=subject sscp rcorr;
	lsmeans Incongruent /adjust=tukey pdiff;
	ods output LSMeans = flankerInteraction;
quit;


proc print data = flankerInteraction;
run;

*interaction plot;
proc plot data = flankerInteraction;
	plot Estimate * Incongruent = CongruentBlock;
quit;

*transpose once again to bring this into one row per person;
proc transpose data = transFlanker out=flankerRTs (rename=(col1=I_congruentMean col2=I_incongruentMean col3 = C_congruentMean col4 = C_incongruentMean));
	var harmean;
	by subject;
run;

proc transpose data = transTrim out=trimFlankerRTs (rename=(col1=T_I_congruentMean col2=T_I_incongruentMean col3 = T_C_congruentMean col4 = T_C_incongruentMean));
	var harmean;
	by subject;
run;

proc print data = flankerRTs (where=(subject=105));
run;

proc print data = trimFlankerRTs (where=(subject=105));
run;

%inc "&projectPath\CorrWithMacro.sas";

proc sort data = diss.questionnaires;
	by subject;
run;

proc sort data = diss.flankerscores;
	by subject;
run;

data flankermerge;
	merge diss.questionnaires diss.flankerscores;
	by subject;
run;

*correlations of RTs;
%corrWith(flankerMerge, I_congruentMean C_congruentMean I_incongruentMean C_incongruentMean, SNAPFactors);
%corrWith(flankerMerge, I_incongruentMean C_incongruentMean, SNAPFacets, partial=C_congruentMean);
%corrWith(flankerMerge, I_congruentMean C_congruentMean I_incongruentMean C_incongruentMean, MPQFactors);
%corrWith(flankerMerge, I_congruentMean C_congruentMean I_incongruentMean C_incongruentMean, MPQSubfactors);
%corrWith(flankerMerge, I_congruentMean C_congruentMean I_incongruentMean C_incongruentMean, MPQFacets);
%corrWith(flankerMerge, I_congruentMean C_congruentMean I_incongruentMean C_incongruentMean, ATQ);

%corrWith(flankerMerge, flankerControl C_flankerControl I_flankerControl, ATQ_ECTotal);
%corrWith(flankerMerge, flankerControl C_flankerControl I_flankerControl, snapfacets, rank=true);


%corrWith(diss.allcombined, I_MeanRT, SNAPFacets, partial=C_MeanRT, rank=true);
%corrWith(diss.allcombined, I_MeanRT_Resid, SNAPFacets);

