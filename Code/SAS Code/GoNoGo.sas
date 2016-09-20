%let projectpath = C:\Documents and Settings\Michael Hallquist\My Documents\Online_Documents\Lab\Dissertation\Analyses;
%let datapath = &projectPath\Data;
*%let datapath = E:\Online_Documents\Lab\Dissertation\Analyses\Data;
libname diss "&datapath";


proc contents data = mergedGoNoGo;
run;

proc print data = mergedGoNoGo;
run;

proc corr data = mergedGoNoGo;
	var MeanGoRT_1 MeanGoRT_3 MeanGoRT_5 MeanGoRT_7
		sumFalseAlarms_1 sumFalseAlarms_3 sumFalseAlarms_5 sumFalseAlarms_7;
run;

proc corr data = diss.allcombined;
	var sumFalseAlarms_1 sumMisses;
run;


proc corr data = mergedGoNoGo;
	var meanRT sumAlarms;
run;

proc univariate data = goNoGo plots;
	var GoDisplay_RT;
run;



*extract just no-go trials;
proc print data = diss.GoNoGoScores;
run;

proc means data = diss.GoNoGoScores;
	class BlockType;
	var falsealarm_sum;
run;

proc mixed data = NumFalseAlarms covtest;
	class BlockType subject;
	model falsealarm_sum = BlockType;
	random intercept /sub=subject;
	lsmeans BlockType /adjust=tukey pdiff;
quit;


ods html;
ods graphics on;
proc glimmix data = noGos plots=studentpanel(unpack);
	class BlockType;
	model falseAlarm (event='1') = BlockType /solution dist=binary link=logit oddsratio;
	random _residual_ /sub=subject type=cs;
	nloptions technique=nrridg;
	lsmeans BlockType /odds oddsratio diff pdiff adjust=tukey plots=meanplot;
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
	lsmestimate BlockType -1 3 -1 -1 /exp e cl divisor=3;  *these estimate ORs > 1, easier to interp;

run;
ods graphics off;
ods html close;

proc means data = noGos std mean stderr;
	class BlockType;
	var falsealarm;
run;

ods html;
ods graphics on;
proc glimmix data = noGos plots=studentpanel(unpack);
	class BlockType;
	model falseAlarm (event='1') = BlockType /solution dist=normal link=identity;
	random _residual_ /sub=subject type=cs;
	nloptions technique=nrridg;
	lsmeans BlockType /diff pdiff adjust=tukey plots=meanplot;
	output out=glimOut pred=predicted residual=resid pearson=pearson student=student;
	lsmestimate BlockType -1 3 -1 -1 /e cl divisor=3;  *these estimate ORs > 1, easier to interp;

run;
ods graphics off;
ods html close;


proc univariate data = NumFalseAlarms plots;
	class BlockType;
	var falsealarm_sum;
run;





proc print data = falseAlarms;
run;

proc sort data = falseAlarms;
	by subject;
run;

proc sort data = diss.questionnaires;
	by subject;
run;

data falseAlarmsMerge;
	merge falseAlarms diss.questionnaires;
	by subject;
run;



data _null_;
	set falseAlarmsMerge;
	file "&dataPath\Mplus\MplusGoNoGo.dat" lrecl=32000;
	if _n_ = 1 then put "subject postemp negtemp disinh mps_per mps_ner mps_cor	onego threego fivego sevengo";
	put subject postemp negtemp disinh mps_per mps_ner mps_cor
		sumFalseAlarms_1 sumFalseAlarms_3 sumFalseAlarms_5 sumFalseAlarms_7;
run;

proc corr data = falseAlarmsMerge;
	var OneGoTrial ThreeGoTrial FiveGoTrial SevenGoTrial;
	with postemp negtemp disinh mps_per mps_ner mps_cor;
run;

proc contents data = diss.questionnaires;
run;


proc means data = goNoGo;
	var GoDisplay_RT;
run;

proc sort data = goNoGo;
	by subject;
run;

proc univariate data = goNoGo plots;
	title "univariate distribution of reaction times";
	title2 "Go-No Go Task, Go Trials";
	title3 "aggregated across block types";
	var GoDisplay_RT;
	by subject;
run;

proc mixed data = goNoGo;
	class BlockType subject;
	model GoDisplay_RT = BlockType;
	random intercept /sub=subject;
	lsmeans BlockType /adjust=tukey pdiff;

quit;

proc print data = numFailedResponses;
run;

proc means data = numFailedResponses;
	class BlockType;
	var failedResponse_Sum;
run;

proc glm data = numFailedResponses;
	class BlockType;
	model failedResponse_Sum = BlockType;
	lsmeans BlockType /adjust=tukey pdiff;
quit;

proc freq data = numFailedResponses;
	tables failedResponse_sum * BlockType;
run;

proc freq data = noGos;
	tables subject * BlockType;
run;

proc print data = noGos (obs=100);
	var subject block NoGoDisplay_OnsetTime NoGoDisplay_RESP falsealarm;
run;

proc print data = noGos (obs=10);
run;

proc freq data = noGos;
	tables subject * falsealarm;
run;



proc print data = NumFalseAlarms;
run;


proc freq data = goNoGo;
	tables BlockType * failedResponse;
run;


proc contents data = diss.GoNoGoScores;
run;



proc univariate data = diss.GoNoGoScores plots;
	var sumFalseAlarms;
	*by subject;
	id subject;
run;

*look for funky RTs;
data findOutliers;
	set diss.GoNoGo;
	if GoDisplay_RT = 0 then GoDisplay_RT = .;
run;

proc sort data = findOutliers;
	by subject;
run;

proc univariate data = findOutliers plots;
	var GoDisplay_RT;
	by subject;
run;

proc print data = diss.GonogoScores;
run;

proc univariate data = diss.GoNoGoScores plots;
	var sumFalseAlarms;
run;

proc sort data = diss.GoNogoScores;
	by subject;
run;

proc sort data = diss.questionnaires;
	by subject;
run;

data merged;
	merge diss.gonogoscores diss.questionnaires;
	by subject;
run;


data _null_;
	set diss.gonogoscores;
	file "&projectPath\Go-No Go\gngscores.dat";
	put subject sumFalsealarms sumMisses sumFalseAlarms_1 sumFalseAlarms_3 sumFalseAlarms_5 sumFalseAlarms_7;
run;

%inc "&projectPath\CorrWithMacro.sas";
*correlations of accuracy scores;
%corrWith(diss.allcombined, MeanGoRT_1 MeanGoRT_3 MeanGoRT_5 MeanGoRT_7, SNAPFacets);
%corrWith(diss.allcombined, sumFalsealarms sumFalseAlarms_1 sumFalseAlarms_3 sumFalseAlarms_5 sumFalseAlarms_7, SNAPFacets, rank=true, partial=sumMisses);
%corrWith(diss.allcombined, sumFalseAlarms_1, SNAPFacets, rank=true, partial=sumMisses);
%corrWith(diss.allcombined, sumFalseAlarms_1, SNAPFacets, rank=true);
%corrWith(diss.allcombined, faresid, SNAPFacets, rank=true);

%corrWith(diss.allcombined, sumFalsealarms sumFalseAlarms_1 sumFalseAlarms_3 sumFalseAlarms_5 sumFalseAlarms_7, SNAPFacets, rank=true, partial=sumMisses);
%corrWith(diss.allcombined, sumFalsealarms sumFalseAlarms_1 sumFalseAlarms_3 sumFalseAlarms_5 sumFalseAlarms_7, SNAPFactors);
%corrWith(diss.allcombined, sumMisses sumMisses_1 sumMisses_3 sumMisses_5 sumMisses_7, SNAPFacets);

%corrWith(diss.allcombined, sumFalsealarms sumFalseAlarms_1, SNAPFacets, rank=true, partial=meanGoRT);

proc corr data = diss.allcombined;
	var sumFalseAlarms_1 sumMisses MeanGoRT;
run;

%corrWith(merged, sumFalsealarms sumFalseAlarms_1 sumFalseAlarms_357, SNAPFacets);

proc corr data = merged ;
	var sumFalseAlarms_1 sumFalseAlarms_357 Inc_PctInaccurate;
run;

proc import datafile="&projectPath\Go-No Go\gngscores.dat" out = work.temp
	dbms=dlm replace;
	getnames=no;
	
run;

data temp;
	set temp;
	drop var1-var5;
	rename var6 = subject
		var7 = gnginhib;
run;

proc sort data = diss.questionnaires;
	by subject;
run;

proc sort data = temp;
	by subject;
run;

data mergeme;
	merge temp diss.questionnaires;
	by subject;
run;

%corrWith(mergeme, gnginhib, SNAPFacets);

proc print data = temp;
	run;

proc print data = temp;
run;


proc import data = temp;
