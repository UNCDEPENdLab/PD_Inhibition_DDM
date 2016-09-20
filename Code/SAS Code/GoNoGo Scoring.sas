*%let projectpath = C:\Documents and Settings\Michael Hallquist\My Documents\Online_Documents\Lab\Dissertation\Analyses;
%let projectpath = C:\Documents and Settings\hallquistmn\My Documents\Online_Documents\Lab\Dissertation\Analyses;
%let datapath = &projectPath\Data;
*%let datapath = E:\Online_Documents\Lab\Dissertation\Analyses\Data;
libname diss "&datapath";

*proc print data = diss.GoNoGo (obs = 100);
*run;

*determine failed responses on go trials;
data goNoGo;
	set diss.GoNoGo;
	if GoDisplay_RT = 0 then do; *zero means no response;
		put "failed response: " subject ", obs: " _n_;
		GoDisplay_RT = .;
		failedResponse = 1;
	end;
	else failedResponse = 0;
	if not missing(GoDisplay_RT) and GoDisplay_RT < 100 then delete; *too fast!;
	if subject = 7 then delete; *person had 32 out of 48 false alarms;
run;

*calculate the sum of failed responses by blocktype;
ods listing close;
proc tabulate data = goNoGo;
	class subject BlockType;
	var failedResponse;
	tables subject * failedResponse * BlockType, sum;
	ods output table=numFailedResponses;
run;
ods listing;

proc sort data = numFailedResponses;
	by subject;
run;

*transpose such that each participant has four scores: one for each blocktype;
proc transpose data = numFailedResponses out = numFailedResponses;
	var failedResponse_sum;
	by subject;
	id BlockType;
run;

*rename variables;
data numFailedResponses;
	set numFailedResponses;
	rename OneGoTrial = sumMisses_1
		ThreeGoTrial = sumMisses_3
		FiveGoTrial = sumMisses_5
		SevenGoTrial = sumMisses_7;
	*rename does not take effect until after data step;
	sumMisses = sum(of OneGoTrial ThreeGoTrial FiveGoTrial SevenGoTrial);
	drop _NAME_;
run;

proc sort data = goNoGo;
	by subject Block;
run;

*CALCULATE NUMBER OF FALSE ALARMS;
data noGos;
	set goNoGo;
	by subject Block;
	if NoGoDisplay_RESP = "{SPACE}" then falsealarm=1;
	else falsealarm=0;
	*no gos are recorded multiple times per go block, so just keep the first in each block;
	if first.Block then output;

run;

proc tabulate data = noGos;
	class subject BlockType;
	tables subject * BlockType;
run;

ods trace on;
*print the sum of the number of false alarms;
ods listing close;
proc tabulate data = noGos;
	class subject BlockType;
	var falsealarm;
	tables subject * falsealarm * BlockType, sum;
	ods output table=falseAlarms;
run;
ods listing;
ods trace off;

proc print data = falseAlarms;
run;

proc mixed data = falseAlarms covtest;
	class BlockType subject;
	model falseAlarm_Sum = BlockType;
	*random intercept /sub=subject;
	repeated /type = un subject=subject sscp rcorr;
	lsmeans BlockType /adjust=tukey pdiff;
	*estimate 'unfamiliar versus familiar' _NAME_ -1 0.333 0.333 0.334 0 /E;
	*estimate 'response interference v. familiar' _NAME_ 0 1 -0.5 -0.5 0 /E;
	estimate 'one versus all others' BlockType -.333 1 -.333 -.334 /E;
quit;

proc mixed data = falseAlarms covtest;
	class BlockType subject;
	model falseAlarm_Sum = BlockType;
	*random intercept /sub=subject;
	repeated /type = cs subject=subject rcorr;
	lsmeans BlockType /adjust=tukey pdiff;
	*estimate 'unfamiliar versus familiar' _NAME_ -1 0.333 0.333 0.334 0 /E;
	*estimate 'response interference v. familiar' _NAME_ 0 1 -0.5 -0.5 0 /E;
	estimate 'one versus all others' BlockType -.333 1 -.333 -.334 /E;
quit;


*collapse false alarm sums into four columns, yielding one row per subject;
proc transpose data = falseAlarms out=falseAlarms;
	var falsealarm_sum;
	by subject;
	id BlockType;
run;

*rename vars and calculate sum;
data falseAlarms;
	set falseAlarms;
	drop _NAME_;
	rename FiveGoTrial = sumFalseAlarms_5
		ThreeGoTrial = sumFalseAlarms_3
		OneGoTrial = sumFalseAlarms_1
		SevenGoTrial = sumFalseAlarms_7;
	sumFalseAlarms = sum(of OneGoTrial ThreeGoTrial FiveGoTrial SevenGoTrial);
run;

*CALCULATE MEAN RT FOR GO TRIALS;
proc univariate data = goNoGo plots noprint;
	class BlockType;
	by subject;
	var GoDisplay_RT;
	output out=IQRCapture qrange=iqr q3=q3; *median=median q1=q1;
run;

/*proc print data = IQRCapture;
run;
*/

proc sort data = IQRCapture;
	by subject;
run;

proc transpose data = IQRCapture out=iqrOut (rename=(col1=iqr_5 col2=iqr_1 col3=iqr_7 col4=iqr_3));
	var iqr;
	by subject;
run;

proc transpose data = IQRCapture out=q3 (rename=(col1=q3_5 col2=q3_1 col3=q3_7 col4=q3_3));
	var q3;
	by subject;
run;

data IQRMerge;
	merge goNogo iqrOut q3;
	by subject;
	drop _NAME_ _LABEL_;
run;

/*proc print data = IQRMerge (obs=500);
	var subject GoDisplay_RT blocktype iqr_1 iqr_3 iqr_5 iqr_7 q3_1 q3_3 q3_5 q3_7;
run;
*/

data GoNoGoTrim;
	set IQRMerge;
	if BlockType = "OneGoTrial" then do;
		if GoDisplay_RT > 3*iqr_1 + q3_1 then GoDisplay_RT = .;
	end;
	else if BlockType = "ThreeGoTrial" then do;
		if GoDisplay_RT > 3*iqr_3 + q3_3 then GoDisplay_RT = .;
	end;
	else if BlockType = "FiveGoTrial" then do;
		if GoDisplay_RT > 3*iqr_5 + q3_5 then GoDisplay_RT = .;
	end;
	else if BlockType = "SevenGoTrial" then do;
		if GoDisplay_RT > 3*iqr_7 + q3_7 then GoDisplay_RT = .;
	end;
run;

proc sort data = goNoGo;
	by subject BlockType;
run;

*flip the reaction times so that they span across columns, stratified by blocktype;
proc transpose data = goNoGo out = goRTs;
	var GoDisplay_RT;
	by subject BlockType;
run;

*proc print data = rtcols;
*run;

*calculate harmonic mean of reaction times;
data goRTs;
	set goRTs;
	harmean = harmean(of col1-col84);
run;

proc sort data = goRTs;
	by subject BlockType;
run;

*collapse blocks into four columns, one row per participant;
proc transpose data = goRTs out = goRTs;
	var harmean;
	by subject;
	id BlockType;
run;

*proc print data = goRTs;
*run;

*rename columns;
data goRTs;
	set goRTs;
	drop _NAME_;
	rename FiveGoTrial = MeanGoRT_5
		ThreeGoTrial = MeanGoRT_3
		SevenGoTrial = MeanGoRT_7
		OneGoTrial = MeanGoRT_1;

	*note that this represents the mean of means, not the mean of all the individual trials;
	MeanGoRT = mean(of OneGoTrial ThreeGoTrial FiveGoTrial SevenGoTrial);

run;


*DO THE SAME FOR TRIMMED;

proc sort data = GoNoGoTrim;
	by subject BlockType;
run;

*flip the reaction times so that they span across columns, stratified by blocktype;
proc transpose data = GoNoGoTrim out = TrimmedGoRTs;
	var GoDisplay_RT;
	by subject BlockType;
run;

*proc print data = rtcols;
*run;

*calculate harmonic mean of reaction times;
data TrimmedGoRTs;
	set TrimmedGoRTs;
	harmean = harmean(of col1-col84);
run;

proc sort data = TrimmedGoRTs;
	by subject BlockType;
run;

*collapse blocks into four columns, one row per participant;
proc transpose data = TrimmedGoRTs out = TrimmedGoRTs;
	var harmean;
	by subject;
	id BlockType;
run;

*proc print data = goRTs;
*run;

*rename columns;
data TrimmedGoRTs;
	set TrimmedGoRTs;
	drop _NAME_;
	rename FiveGoTrial = T_MeanGoRT_5
		ThreeGoTrial = T_MeanGoRT_3
		SevenGoTrial = T_MeanGoRT_7
		OneGoTrial = T_MeanGoRT_1;

	*note that this represents the mean of means, not the mean of all the individual trials;
	T_MeanGoRT = mean(of OneGoTrial ThreeGoTrial FiveGoTrial SevenGoTrial);

run;



*merge calculated scores for Go-No Go;
proc sort data = goRTs;
	by subject;
run;

proc sort data = TrimmedGoRTs;
	by subject;
run;

proc sort data = falseAlarms;
	by subject;
run;

proc sort data = numFailedResponses;
	by subject;
run;

data diss.GoNoGoScores;
	merge TrimmedGoRTs goRTs falseAlarms numFailedResponses;
	by subject;
	if subject = 7 then delete; *person had 32 out of 48 false alarms;
run;


proc reg data = diss.GoNoGoScores;
	model sumFalseAlarms_1 = sumMisses;
	output out=faresid residual=faresid;
quit;

proc sort data = faresid;
	by subject;
run;

proc sort data = diss.GoNoGoScores;
	by subject;
run;

data diss.GoNoGoScores;
	merge diss.GoNoGoScores faresid (keep=subject faresid);
	by subject;
run;

proc print data = diss.goNogoscores;
run;
/*
proc means data = diss.goNogoscores;
	var MeanGoRT_1 T_MeanGoRT_1 MeanGoRT_3 T_MeanGoRT_3 MeanGoRT_5 T_MeanGoRT_5 MeanGoRT_7 T_MeanGoRT_7 ;
run;
*/
