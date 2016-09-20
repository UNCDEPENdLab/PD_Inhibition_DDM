%let projectpath = C:\Documents and Settings\Michael Hallquist\My Documents\Online_Documents\Lab\Dissertation\Analyses;
%let datapath = &projectPath\Data;
*%let datapath = E:\Online_Documents\Lab\Dissertation\Analyses\Data;
libname diss "&datapath";
options fmtsearch=(diss);

%MACRO EmotFlankRTReliability();
	proc sort data = diss.emotFlanker;
		by subject colidentifier;
	run;

	proc transpose data = diss.emotFlanker out=emotReliability; *(rename=(col1=I_congruentMean col2=C_congruentMean col3 = I_incongruentMean col4 = C_incongruentMean));
		var FlankerDisplay_RT;
		by subject colidentifier;
	run;

	proc sort data = emotReliability;
		by colidentifier;
	run;

	proc corr data = emotReliability nocorr alpha;
		var col1-col12;
		by colidentifier;
	run;
%MEND;



*combine means for congruent trials for the same target emotion across blocks;
proc means data = diss.EmotFlankerRTScores;
	var C_AngHap_Ang C_AngNeu_Ang C_AngSad_Ang;
run;

proc glm data = diss.EmotFlankerRTScores;
	model C_AngHap_Ang C_AngNeu_Ang C_AngSad_Ang = /nouni;
	repeated condition 3 profile /PRINTE;
quit;

proc ttest data = diss.EmotFlankerRTScores;
	paired C_AngHap_Ang*C_AngNeu_Ang C_AngNeu_Ang*C_AngSad_Ang;
run;

proc factor data=diss.EmotFlankerRTScores 
	priors=smc msa scree residual preplot
	rotate=promax reorder plot 
	outstat=fact_all;
	title3 'Principal Factor Analysis with Promax Rotation';
	var I_AngHap_Ang I_AngHap_Hap I_AngNeu_Ang I_AngNeu_Neu I_AngSad_Sad 
		I_AngSad_Ang I_HapNeu_Hap I_HapNeu_Neu I_SadHap_Sad I_SadHap_Hap I_SadNeu_Sad I_SadNeu_Neu
		C_AngHap_Ang C_AngHap_Hap C_AngNeu_Ang C_AngNeu_Neu C_AngSad_Sad 
		C_AngSad_Ang C_HapNeu_Hap C_HapNeu_Neu C_SadHap_Sad C_SadHap_Hap C_SadNeu_Sad C_SadNeu_Neu;
run;


proc factor data=diss.EmotFlankerRTScores simple corr;
	title3 'Principal Component Analysis';
	var I_AngHap_Ang I_AngHap_Hap I_AngNeu_Ang I_AngNeu_Neu I_AngSad_Sad 
	I_AngSad_Ang I_HapNeu_Hap I_HapNeu_Neu I_SadHap_Sad I_SadHap_Hap I_SadNeu_Sad I_SadNeu_Neu
	C_AngHap_Ang C_AngHap_Hap C_AngNeu_Ang C_AngNeu_Neu C_AngSad_Sad 
	C_AngSad_Ang C_HapNeu_Hap C_HapNeu_Neu C_SadHap_Sad C_SadHap_Hap C_SadNeu_Sad C_SadNeu_Neu;
run;


proc glm data = diss.EmotFlankerRTScores;
	model I_AngHap_Ang I_AngNeu_Ang I_AngSad_Ang = ;
	repeated condition 3 /PRINTE;
	lsmeans condition;
quit;


proc print data = diss.EmotFlankerRTScores (where=(subject=114));
run;

proc corr data = diss.EmotFlankerRTScores rank;
	var AngHap_Ang_Diff AngHap_Hap_Diff AngNeu_Ang_Diff AngNeu_Neu_Diff AngSad_Sad_Diff
		AngSad_Ang_Diff HapNeu_Hap_Diff HapNeu_Neu_Diff SadHap_Sad_Diff SadHap_Hap_Diff
		SadNeu_Sad_Diff SadNeu_Neu_Diff;

run;

proc corr data = diss.EmotFlankerRTScores rank;
	var I_AngHap_Ang I_AngHap_Hap I_AngNeu_Ang I_AngNeu_Neu I_AngSad_Sad 
		I_AngSad_Ang I_HapNeu_Hap I_HapNeu_Neu I_SadHap_Sad I_SadHap_Hap I_SadNeu_Sad I_SadNeu_Neu;
run;


ods trace on;
proc corr data = diss.EmotFlankerRTScores; *rank;
	var C_AngHap_Ang C_AngHap_Hap C_AngNeu_Ang C_AngNeu_Neu C_AngSad_Sad 
		C_AngSad_Ang C_HapNeu_Hap C_HapNeu_Neu C_SadHap_Sad C_SadHap_Hap C_SadNeu_Sad C_SadNeu_Neu;
	ods output PearsonCorr=matrixOut;
run;
ods trace off;

proc print data = matrixout;
run;

proc contents data = matrixout;
run;

data matrixOut;
	set matrixOut;
	drop PC_AngHap_Ang PC_AngHap_Hap PC_AngNeu_Ang PC_AngNeu_Neu PC_AngSad_Ang
		PC_AngSad_Sad PC_HapNeu_Hap PC_HapNeu_Neu PC_SadHap_Hap PC_SadHap_Sad 
		PC_SadNeu_Neu PC_SadNeu_Sad;
run;

proc sort data = matrixOut;
	by variable;
run;

proc transpose data =matrixOut out=matrixOut;
	by variable;
run;

proc print data = matrixOut;
	where col1 < 1.0 and col1 > 0.7;
run;

*the results from the correlations above suggest that RT for congruent trials is highly
	correlated within the block. For example, for the sad-anger block, RTs for congruent
	sad and angry trials are highly correlated. Thus, look at congruent RTs within blocks;

proc corr data = diss.EmotFlankerRTScores;
	var C_AngSad_Sad C_AngSad_Ang;
run;

proc corr data = diss.EmotFlankerRTScores;
	var C_HapNeu_Hap C_HapNeu_Neu;
run;

proc corr data = diss.EmotFlankerRTScores;
	var C_SadNeu_Sad C_SadNeu_Neu;
run;

proc corr data = diss.EmotFlankerRTScores;
	var C_SadHap_Sad C_SadHap_Hap;
run;

proc corr data = diss.EmotFlankerRTScores;
	var C_AngNeu_Ang C_AngNeu_Neu;
run;

proc corr data = diss.EmotFlankerRTScores;
	var C_AngHap_Ang C_AngHap_Hap;
run;

*add incongruent rts in the mix;
proc corr data = diss.EmotFlankerRTScores;
	var C_AngSad_Sad C_AngSad_Ang I_AngSad_Sad I_AngSad_Ang;
run;

proc corr data = diss.EmotFlankerRTScores;
	var C_HapNeu_Hap C_HapNeu_Neu I_HapNeu_Hap I_HapNeu_Neu;
run;

proc corr data = diss.EmotFlankerRTScores;
	var C_SadNeu_Sad C_SadNeu_Neu I_SadNeu_Sad I_SadNeu_Neu;
run;

proc corr data = diss.EmotFlankerRTScores;
	var C_SadHap_Sad C_SadHap_Hap I_SadHap_Sad I_SadHap_Hap;
run;

proc corr data = diss.EmotFlankerRTScores;
	var C_AngNeu_Ang C_AngNeu_Neu I_AngNeu_Ang I_AngNeu_Neu;
run;

proc corr data = diss.EmotFlankerRTScores;
	var C_AngHap_Ang C_AngHap_Hap I_AngNeu_Ang I_AngNeu_Neu;
run;

proc sort data = diss.EmotFlankerRTScores;
	by subject;
run;

proc sort data = diss.questionnaires;
	by subject;
run;

proc sort data = diss.EmotFlankerACCScores;
	by subject;
run;

data emotFlankMerge;
	merge diss.EmotFlankerRTScores diss.EmotFlankerACCScores diss.questionnaires;
	by subject;
run;

proc corr data = emotFlankMerge rank;
	var AngHap_Ang_Diff AngHap_Hap_Diff AngNeu_Ang_Diff AngNeu_Neu_Diff;
	with NEGTEMP MISTRUST MANIP AGG SELFHARM ECCPERC DEPEN POSTEMP 
	EXHIB ENTITL DETACH DISINH IMPUL PROPER HDWK LOSLFEST SUICPRON
	MPS_PER MPS_NER MPS_COR;
run;

data _null_;
	set emotFlankMerge;
	file "&dataPath\Mplus\MplusEmotFlanker.dat" lrecl=32000;
	if _n_ = 1 then put "subject postemp negtemp disinh mps_per mps_ner mps_cor	I_AngHap_Ang I_AngHap_Hap I_AngNeu_Ang I_AngNeu_Neu I_AngSad_Sad I_AngSad_Ang I_HapNeu_Hap I_HapNeu_Neu I_SadHap_Sad I_SadHap_Hap I_SadNeu_Sad I_SadNeu_Neu C_AngHap_Ang C_AngHap_Hap C_AngNeu_Ang C_AngNeu_Neu C_AngSad_Sad C_AngSad_Ang C_HapNeu_Hap C_HapNeu_Neu C_SadHap_Sad C_SadHap_Hap C_SadNeu_Sad C_SadNeu_Neu";
	put subject postemp negtemp disinh mps_per mps_ner mps_cor
		I_AngHap_Ang I_AngHap_Hap I_AngNeu_Ang I_AngNeu_Neu I_AngSad_Sad 
		I_AngSad_Ang I_HapNeu_Hap I_HapNeu_Neu I_SadHap_Sad I_SadHap_Hap 
		I_SadNeu_Sad I_SadNeu_Neu C_AngHap_Ang C_AngHap_Hap C_AngNeu_Ang C_AngNeu_Neu C_AngSad_Sad 
		C_AngSad_Ang C_HapNeu_Hap C_HapNeu_Neu C_SadHap_Sad C_SadHap_Hap C_SadNeu_Sad C_SadNeu_Neu;
run;


%inc "&projectPath\CorrWithMacro.sas";

options mprint;
%corrWith(diss.allcombined, I_PctInaccurate C_PctInaccurate Ang_PctInaccurate Hap_PctInaccurate Neu_PctInaccurate Sad_PctInaccurate PctInaccurate, SNAPFactors, rank=true);

%corrWith(diss.allcombined, I_PctInaccurate C_PctInaccurate Ang_PctInaccurate Hap_PctInaccurate Neu_PctInaccurate Sad_PctInaccurate PctInaccurate, MPQFactors);
%corrWith(diss.allcombined, I_PctInaccurate C_PctInaccurate Ang_PctInaccurate Hap_PctInaccurate Neu_PctInaccurate Sad_PctInaccurate PctInaccurate, MPQFacets);
%corrWith(diss.allcombined, I_PctInaccurate C_PctInaccurate Ang_PctInaccurate Hap_PctInaccurate Neu_PctInaccurate Sad_PctInaccurate PctInaccurate, ATQ);

%corrWith(diss.allcombined, I_PctInaccurate, SNAPFacets, rank=true, partial=C_PctInaccurate);

%corrWith(diss.allcombined, I_PctInaccurate, SNAPFacets, rank=true);

ods html;
%corrWith(diss.allcombined, I_PctInaccurate C_PctInaccurate Ang_PctInaccurate Hap_PctInaccurate Neu_PctInaccurate Sad_PctInaccurate PctInaccurate, SNAPFacets, rank=true);
%corrWith(diss.allcombined, I_PctInaccurate, SNAPFacets, rank=true, partial=C_PctInaccurate);
%corrWith(diss.allcombined, Ang_I_PctInaccurate Sad_I_PctInaccurate Ang_C_PctInaccurate Sad_C_PctInaccurate, SNAPFacets, rank=true);
%corrWith(diss.allcombined, AngSad_I_Sad_PctInaccurate AngSad_I_Ang_PctInaccurate, SNAPFacets, rank=true);

*THESE 5 ARE THE KEEPERS (SEE EVERNOTE);
%corrWith(diss.allcombined, Sad_I_PctInaccurate, SNAPFacets, rank=true, partial=Sad_C_PctInaccurate);
%corrWith(diss.allcombined, Ang_I_PctInaccurate, SNAPFacets, rank=true, partial=Ang_C_PctInaccurate);
%corrWith(diss.allcombined, I_PctInaccurate, SNAPFacets, rank=true, partial=C_PctInaccurate);
%corrWith(diss.allcombined, I_Ang, SNAPFacets, rank=true, partial=C_Ang);
%corrWith(diss.allcombined, I_Sad, SNAPFacets, rank=true, partial=C_Sad);

proc corr data = diss.allcombined;
	var Sad_I_PctInaccurate Sad_C_PctInaccurate;
run;

proc corr data = diss.allcombined;
	var Ang_I_PctInaccurate Ang_C_PctInaccurate;
run;

proc corr data = diss.allcombined;
	var I_PctInaccurate C_PctInaccurate;
run;

proc corr data = diss.allcombined;
	var I_Sad C_Sad;
run;

proc corr data = diss.allcombined;
	var I_Ang C_Ang;
run;


%corrWith(diss.allcombined, Hap_I_PctInaccurate, SNAPFacets, rank=true, partial=Hap_C_PctInaccurate);
%corrWith(diss.allcombined, Neu_I_PctInaccurate, SNAPFacets, rank=true, partial=Neu_C_PctInaccurate);

%corrWith(diss.allcombined, Sad_I_PctInaccurate, SNAPFacets, rank=true, partial=C_PctInaccurate);
%corrWith(diss.allcombined, Ang_I_PctInaccurate, SNAPFacets, rank=true, partial=C_PctInaccurate);
%corrWith(diss.allcombined, Hap_I_PctInaccurate, SNAPFacets, rank=true, partial=C_PctInaccurate);
%corrWith(diss.allcombined, Neu_I_PctInaccurate, SNAPFacets, rank=true, partial=C_PctInaccurate);

%corrWith(diss.allcombined, AngSad_I_Ang_PctInaccurate, SNAPFacets, rank=true, partial=AngSad_C_Ang_PctInaccurate);
%corrWith(diss.allcombined, AngSad_I_Sad_PctInaccurate, SNAPFacets, rank=true, partial=AngSad_C_Sad_PctInaccurate);

proc contents data = diss.EmotFlankerRTScores;
run;
ods html close;

proc corr data = diss.allcombined;
	var I_PctInaccurate C_PctInaccurate Ang_PctInaccurate Hap_PctInaccurate Neu_PctInaccurate Sad_PctInaccurate PctInaccurate;
run;


proc corr data = diss.allcombined;
	var Ang_I_PctInaccurate Ang_C_PctInaccurate Sad_I_PctInaccurate Sad_C_PctInaccurate;
run;

proc corr data = diss.allcombined;
	var AngSad_I_Sad_PctInaccurate AngSad_I_Ang_PctInaccurate AngSad_C_Sad_PctInaccurate AngSad_C_Ang_PctInaccurate;
run;

proc univariate data = diss.allcombined;
	var Ang_C_PctInaccurate Sad_C_PctInaccurate Hap_C_PctInaccurate C_PctInaccurate;
run;

proc contents data = diss.emotFlanker;
run;


proc univariate data = diss.emotFlanker plots;
	var FlankerDisplay_RT;
run;





proc freq data = diss.emotFlanker;
	tables condition;
run;

proc freq data = diss.emotFlanker;
	tables emotion * congruent;
run;

proc contents data = diss.emotFlanker;
run;

proc means data = diss.emotFlanker;
	class condition;
	var FlankerDisplay_RT;
run;

proc sort data = diss.emotFlanker;
	by subject;
run;

proc freq data = diss.emotFlanker;
	tables FlankerDisplay_ACC * emotion * condition;
	by subject;
run;


proc means data = diss.emotFlanker;
	class condition emotion congruent;
	var FlankerDisplay_RT;
run;

data temp;
	set diss.emotFlanker;
	if subject in (98, 100, 40) then delete;
run;

ods html;
ods graphics on;
proc mixed data = temp covtest;
	title 'RT Analysis of Emotional Flanker by emotion, pair, and congruence';
	class emotion subject condition congruent;
	model FlankerDisplay_RT = emotion | congruent | condition(emotion);
	repeated /sub=subject type=cs rcorr group=emotion;
	*test the effect of condition for each emotion;
	lsmeans condition(emotion) /diff adjust=tukey slice=emotion;
	lsmeans emotion /diff adjust=tukey;
	lsmeans congruent /diff adjust=tukey;
	*output out=glimOut pred=predicted residual=resid pearson=pearson student=student;

	format emotion Emotions. condition EmotionPair.;
quit;
ods graphics off;
ods html close;

*FINAL 3-WAY RT MODEL!;
ods html;
ods graphics on;
proc glimmix data = temp plots=studentpanel(unpack) ic=pq;
	class condition emotion congruent subject;
	model FlankerDisplay_RT = condition(emotion) | emotion | congruent /solution dist=normal link=identity;
	random _residual_ /sub=subject type=cs;
	nloptions technique=nrridg;
	lsmeans condition(emotion) /diff slice=emotion slicediff=emotion adjust=tukey plots=meanplot(sliceby=emotion join);
	lsmeans emotion /diff adjust=tukey plots=meanplot(sliceby=emotion join);
	lsmeans congruent /diff adjust=tukey plots=meanplot(sliceby=congruent join);
	lsmeans emotion*congruent /plots=meanplot(sliceby=congruent join);
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
	format emotion Emotions.;

run;
ods graphics off;
ods html close;

*CALCULATE DESCRIPTIVE STATS FOR EMOTIONS;

proc means data = temp;
	class congruent;
	var FlankerDisplay_RT;
run;

proc sort data = temp;
	by subject emotion;
run;

proc transpose data = temp out=byEmotion;
	var FlankerDisplay_RT;
	by subject emotion;
run;

*proc print data = byEmotion;
*run;

data byEmotion;
	set byEmotion;
	harmean = harmean(of col1-col72);
	arithmean = mean (of col1-col72);
run;

proc transpose data = byEmotion out=byEmotion;
	var harmean;
	by subject;
	id emotion;
run;

proc means data = byEmotion mean std;
	var Sad Happy Neutral Angry;
run;



ods rtf file = "&projectPath\Mixed Analysis by Emotion.rtf";
proc mixed data = diss.emotFlanker covtest;
	title 'mean reaction time differences across conditions';
	class emotion subject;
	model FlankerDisplay_RT = emotion;
	repeated /sub=subject type=cs rcorr group=emotion;
	lsmeans emotion;
	format emotion Emotions.;
quit;
ods rtf close;


proc mixed data = diss.emotFlanker (where=(emotion=1)) covtest;
	title 'mean reaction time differences for Sadness';
	class condition congruent subject;
	model FlankerDisplay_RT = condition | congruent;
	repeated /sub=subject type=cs rcorr group=condition*congruent;
	lsmeans condition * congruent;
	format condition EmotionPair.;
quit;


proc mixed data = diss.emotFlanker (where=(emotion=1)) covtest;
	title 'mean reaction time differences for Sadness';
	class condition congruent subject;
	model FlankerDisplay_RT = condition | congruent;
	repeated /sub=subject type=cs rcorr group=condition;
	lsmeans condition * congruent;
	format condition EmotionPair.;
quit;


proc mixed data = diss.emotFlanker (where=(emotion=1)) covtest;
	title 'mean reaction time differences for Sadness';
	class condition congruent subject;
	model FlankerDisplay_RT = condition | congruent;
	repeated /sub=subject type=cs rcorr;
	lsmeans condition * congruent;
	format condition EmotionPair.;
quit;


proc mixed data = diss.emotFlanker (where=(emotion=1)) covtest;
	title 'mean reaction time differences for Sadness';
	class condition congruent subject;
	model FlankerDisplay_RT = condition | congruent;
	repeated /sub=subject type=cs r;
	lsmeans condition * congruent;
	format condition EmotionPair.;
quit;


proc mixed data = diss.emotFlanker (where=(emotion=1)) covtest;
	title 'mean reaction time differences for Sadness';
	class condition subject;
	model FlankerDisplay_RT = condition;
	random intercept /sub=subject type=cs g gcorr;
	format condition EmotionPair.;
quit;




%MACRO RT_ByConditionCongruent(dsn=diss.emotFlanker, emotion=, emotionLabel=);

	proc mixed data = diss.emotFlanker (where=(emotion=&emotion));
		title "mean reaction time differences for &emotionLabel";
		class condition congruent subject;
		model FlankerDisplay_RT = condition | congruent;
		repeated /sub=subject;
		lsmeans condition * congruent;
		lsmeans condition /pdiff adjust=tukey;
		format condition EmotionPair.;
	quit;

%MEND;

%RT_ByConditionCongruent(emotion=1, emotionLabel=Sadness);
%RT_ByConditionCongruent(emotion=2, emotionLabel=Anger);
%RT_ByConditionCongruent(emotion=3, emotionLabel=Happiness);
%RT_ByConditionCongruent(emotion=4, emotionLabel=Neutral);


ods html;
ods graphics on;
proc glimmix data = temp plots=studentpanel(unpack);
	title 'two-way GLMM testing effect of emotion and congruence';
	class emotion congruent;
	model FlankerDisplay_ACC  = emotion | congruent /solution dist=binary link=logit oddsratio;
	random _residual_ /sub=subject type=cs;
	nloptions technique=nrridg;
	lsmeans emotion /odds oddsratio diff adjust=tukey plots=meanplot(sliceby=congruent join);
	lsmeans congruent /odds oddsratio diff adjust=tukey plots=meanplot(sliceby=congruent join);
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
	lsmestimate congruent -1 1 /exp cl;  *these estimate ORs > 1, easier to interp;
	format emotion Emotions.;

run;
ods graphics off;
ods html close;

ods html;
ods graphics on;
proc glimmix data = temp plots=studentpanel(unpack);
	class condition emotion;
	model FlankerDisplay_ACC = condition emotion(condition) /solution dist=binary link=logit oddsratio;
	random _residual_ /sub=subject type=cs;
	nloptions technique=nrridg;
	*lsmeans emotion(condition) /odds oddsratio diff plots=meanplot(sliceby=congruent join);
	*lsmeans congruent /odds oddsratio diff plots=meanplot(sliceby=congruent join);
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
	*lsmestimate congruent -1 1 /exp cl;  *these estimate ORs > 1, easier to interp;
	format emotion Emotions.;

run;
ods graphics off;
ods html close;

*FINAL 3-WAY MODEL WITH BLOCKED COVARIANCE PER EMOTION;
ods html;
ods graphics on;
proc glimmix data = temp plots=studentpanel(unpack);
	class condition emotion congruent subject;
	model FlankerDisplay_ACC = condition(emotion) | emotion | congruent /solution dist=binary link=logit oddsratio;
	random _residual_ /sub=subject type=cs group=emotion;
	nloptions technique=nrridg;
	lsmeans condition(emotion) /odds oddsratio diff slice=emotion slicediff=emotion plots=meanplot(sliceby=emotion join);
	lsmeans emotion /odds oddsratio diff adjust=tukey plots=meanplot(sliceby=emotion join);
	lsmeans congruent /odds oddsratio diff adjust=tukey plots=meanplot(sliceby=congruent join);
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
	*lsmestimate congruent -1 1 /exp cl;  *these estimate ORs > 1, easier to interp;
	format emotion Emotions.;

run;
ods graphics off;
ods html close;

*FINAL 3-WAY MODEL W/O BLOCKING;
ods html;
ods graphics on;
proc glimmix data = temp plots=studentpanel(unpack) ic=pq;
	class condition emotion congruent subject;
	model FlankerDisplay_ACC = condition(emotion) | emotion | congruent /solution dist=binary link=logit oddsratio;
	random _residual_ /sub=subject type=cs;
	nloptions technique=nrridg;
	lsmeans condition(emotion) /odds oddsratio diff slice=emotion slicediff=emotion adjust=tukey plots=meanplot(sliceby=emotion join);
	lsmeans emotion /odds oddsratio diff adjust=tukey plots=meanplot(sliceby=emotion join);
	lsmeans congruent /odds oddsratio diff adjust=tukey plots=meanplot(sliceby=congruent join);
	lsmeans emotion*congruent /plots=meanplot(sliceby=congruent join);
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
	*lsmestimate congruent -1 1 /exp cl;  *these estimate ORs > 1, easier to interp;
	format emotion Emotions.;

run;
ods graphics off;
ods html close;



ods html;
ods graphics on;
proc glimmix data = temp plots=studentpanel(unpack) ic=pq;
	class condition emotion congruent subject;
	model FlankerDisplay_ACC = condition(emotion) | emotion | congruent /solution dist=norma link=identity;
	random _residual_ /sub=subject type=cs;
	nloptions technique=nrridg;
	lsmeans condition(emotion) /diff slice=emotion slicediff=emotion adjust=tukey plots=meanplot(sliceby=emotion join);
	lsmeans emotion /diff adjust=tukey plots=meanplot(sliceby=emotion join);
	lsmeans congruent /diff adjust=tukey plots=meanplot(sliceby=congruent join);
	lsmeans emotion*congruent /plots=meanplot(sliceby=congruent join);
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
	*lsmestimate congruent -1 1 /exp cl;  *these estimate ORs > 1, easier to interp;
	format emotion Emotions.;

run;
ods graphics off;
ods html close;




*residuals analysis to look for problematic outliers;
proc sort data = glimOut;
	by subject;
run;

data residCheck;
	set glimOut;
	where student > 3;
run;

proc freq data = residCheck;
	tables subject /out=outFreq;
run;
proc sort data = outFreq;
	by count;
run;
proc print data = outFreq;
	var subject count;
run;

proc sort data = accuracy;
	by pctinaccurate;
run;

proc print data = accuracy;
	var subject pctinaccurate;
run;

proc freq data = glimOut;
	tables subject;
run;


proc print data = glimOut;
	var subject student resid predicted flankerdisplay_acc;
	where student > 4;
run;

*LOOK AT INACCURACY RATES ACROSS EMOTIONS;
*NOTE THAT THE PER-EMOTION ANALYSES ARE DEPRECATED, AS THE 3-WAY GLMM SUBSUMES THESE;

proc freq data = temp (where=(emotion=1));
	tables condition * congruent * FlankerDisplay_ACC;
	tables condition * FlankerDisplay_ACC;
	tables congruent * FlankerDisplay_ACC;
run;

*analyze differences in sad inaccuracies across blocks;
ods html;
ods graphics on;
proc glimmix data = temp (where=(emotion=1)) plots=studentpanel(unpack);
	class condition congruent;
	model FlankerDisplay_ACC  = condition | congruent /solution dist=binary link=logit oddsratio;
	random _residual_ /sub=subject type=cs gcorr vcorr group=condition;
	nloptions technique=nrridg;
	lsmeans condition /odds oddsratio diff adjust=tukey plots=meanplot(sliceby=congruent join);
	lsmeans congruent /odds oddsratio diff adjust=tukey plots=meanplot(sliceby=congruent join);
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
	lsmestimate congruent -1 1 /exp cl;  *these estimate ORs > 1, easier to interp;
	format condition EmotionPair.;

run;
ods graphics off;
ods html close;


*analyze differences in angry inaccuracies across blocks;
proc freq data = temp (where=(emotion=2));
	tables condition * congruent * FlankerDisplay_ACC;
	tables condition * FlankerDisplay_ACC;
	tables congruent * FlankerDisplay_ACC;
run;


ods html;
ods graphics on;
proc glimmix data = temp (where=(emotion=2)) plots=studentpanel(unpack);
	class condition congruent;
	model FlankerDisplay_ACC  = condition | congruent /solution dist=binary link=logit oddsratio;
	random _residual_ /sub=subject type=cs gcorr vcorr group=condition;
	nloptions technique=nrridg;
	lsmeans condition /odds oddsratio diff adjust=tukey plots=meanplot(sliceby=congruent join);
	lsmeans congruent /odds oddsratio diff adjust=tukey plots=meanplot(sliceby=congruent join);
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
	lsmestimate congruent 'incon v. con' -1 1 /exp cl;  *these estimate ORs > 1, easier to interp;
	format condition EmotionPair.;

run;
ods graphics off;
ods html close;



*analyze differences in happy inaccuracies across blocks;
proc freq data = temp (where=(emotion=3));
	title 'inaccurate responses for happy trials as a function of congruence and emotion pair';
	tables condition * congruent * FlankerDisplay_ACC;
	tables condition * FlankerDisplay_ACC;
	tables congruent * FlankerDisplay_ACC;
run;


ods html;
ods graphics on;
proc glimmix data = temp (where=(emotion=3)) plots=studentpanel(unpack);
	class condition congruent;
	model FlankerDisplay_ACC  = condition | congruent /solution dist=binary link=logit oddsratio;
	random _residual_ /sub=subject type=cs gcorr vcorr group=condition;
	nloptions technique=nrridg;
	lsmeans condition /odds oddsratio diff adjust=tukey plots=meanplot(sliceby=congruent join);
	lsmeans congruent /odds oddsratio diff adjust=tukey plots=meanplot(sliceby=congruent join);
	lsmeans congruent*condition /odds oddsratio diff slice=condition adjust=tukey plots=meanplot(sliceby=congruent join);
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
	lsmestimate congruent 'incon v. con' -1 1 /exp cl;  *these estimate ORs > 1, easier to interp;
	format condition EmotionPair.;

run;
ods graphics off;
ods html close;


*analyze differences in neutral inaccuracies across blocks;
proc freq data = temp (where=(emotion=4));
	title 'inaccurate responses for neutral trials as a function of congruence and emotion pair';
	tables condition * congruent * FlankerDisplay_ACC;
	tables condition * FlankerDisplay_ACC;
	tables congruent * FlankerDisplay_ACC;
run;


ods html;
ods graphics on;
proc glimmix data = temp (where=(emotion=4)) plots=studentpanel(unpack);
	class condition congruent;
	model FlankerDisplay_ACC  = condition | congruent /solution dist=binary link=logit oddsratio;
	random _residual_ /sub=subject type=cs gcorr vcorr group=condition;
	nloptions technique=nrridg;
	lsmeans condition /odds oddsratio diff adjust=tukey plots=meanplot(sliceby=congruent join);
	lsmeans congruent /odds oddsratio diff adjust=tukey plots=meanplot(sliceby=congruent join);
	lsmeans congruent*condition /odds oddsratio diff slice=condition adjust=tukey plots=meanplot(sliceby=congruent join);
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
	lsmestimate congruent 'incon v. con' -1 1 /exp cl;  *these estimate ORs > 1, easier to interp;
	format condition EmotionPair.;

run;
ods graphics off;
ods html close;

