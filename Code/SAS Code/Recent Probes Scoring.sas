%let projectpath = C:\Documents and Settings\Michael Hallquist\My Documents\Online_Documents\Lab\Dissertation\Analyses;
*%let projectpath = E:\Online_Documents\Lab\Dissertation\Analyses;
%let datapath = &projectPath\Data;
libname diss "&datapath";

*proc contents data = diss.RecentProbes;
*run;

*proc print data = diss.RecentProbes (obs=100);
*	var Running_Block_ Running_Trial_;
*run;

*proc freq data = diss.RecentProbes;
*	tables condition;
*run;

*proc univariate data = diss.RecentProbes;
*	var ProbeDisplay_RT;
*run;


data recentProbes;
	set diss.RecentProbes;
	if ProbeDisplay_ACC = 0 then ProbeDisplay_RT = .;
run;

proc sort data = recentProbes;
	by subject condition;
run;


proc transpose data = recentProbes out=probesTrans;
	var ProbeDisplay_RT;
	by subject condition;
run;

proc sort data = probesTrans;
	by condition;
run;

*proc print data = probesTrans;
*run;

proc corr data = probesTrans alpha nocorr;
	title "cronbach alpha values for each condition";
	var Col1-Col9;
	by condition;
run;

data probesTrans;
	set probesTrans;
	harmean = harmean (of col1-col36); *note that only Y trials have 36 trials (the rest have 9);
	arithmean = mean (of col1-col36); *note that only Y trials have 36 trials (the rest have 9);
run;

*proc print data = probesTrans (obs=100);
*run;

proc sort data = probesTrans;
	by subject;
run;


proc transpose data = probesTrans out=arith;
	var arithmean;
	by subject;
	id condition;
run;

proc means data = arith (where=(subject not in (40, 65, 100))); *near 50% in terms of accuracy;;
	title 'arithmetic mean reaction times by condition';
	var Y Nfam0 Nfam1 Nfam1RI Nfam2;
run;

proc transpose data = probesTrans out=probesTrans;
	var harmean;
	by subject;
	id condition;
run;


*proc print data = probesTrans;
*run;

data probesTrans;
	set probesTrans;
	drop _NAME_;
	rename Y = RP_Y_RT
		Nfam0 = RP_Nfam0_RT
		Nfam1 = RP_Nfam1_RT
		Nfam1RI = RP_Nfam1RI_RT
		Nfam2 = RP_Nfam2_RT;

	RP_N_RT = mean (of Nfam0 Nfam1 Nfam1RI Nfam2); *this is a mean of means, not agg of the trials;
	RP_familiar_RT = mean (of Nfam1 Nfam2); *Nfam1RI this is a mean of means, not agg of the trials;
	RP_RT = mean (of RP_N_RT Y);
run;

*CONSIDER SHIFTING SCORING TO SUM, NOT MEAN;
ods listing close;
proc tabulate data = diss.RecentProbes;
	class subject condition;
	var ProbeDisplay_ACC;
	tables subject * ProbeDisplay_ACC * condition, mean;
	ods output table=propAccurate;
run;
ods listing;

*proc print data = propAccurate;
*run;

proc sort data = propAccurate;
	by subject;
run;

proc transpose data = propAccurate out = propAccurate;
	var ProbeDisplay_ACC_MEAN;
	by subject;
	id condition;
run;

*proc print data = propAccurate;
*run;

data propAccurate;
	set propAccurate;
	drop _NAME_;
	rename Y = RP_Y_ACC
		Nfam0 = RP_Nfam0_ACC
		Nfam1 = RP_Nfam1_ACC
		Nfam1RI = RP_Nfam1RI_ACC
		Nfam2 = RP_Nfam2_ACC;
	RP_N_ACC = mean(of Nfam0 Nfam1 Nfam1RI Nfam2);
	RP_ACC = mean(of Y RP_N_ACC);
run;

*proc print data = propAccurate;
*run;

*score error sums;
ods listing close;
proc tabulate data = diss.RecentProbes;
	class subject condition;
	var ProbeDisplay_ACC;
	tables subject * ProbeDisplay_ACC * condition, sum;
	ods output table=sumAccurate;
run;
ods listing;

proc sort data = sumAccurate;
	by subject;
run;

proc transpose data = sumAccurate out = sumAccurate;
	var ProbeDisplay_ACC_SUM;
	by subject;
	id condition;
run;

data sumAccurate;
	set sumAccurate;
	drop _NAME_;
	*reverse the scoring to score for number of errors, not hits;
	Y = 36 - Y;
	Nfam0 = 9 - Nfam0;
	Nfam1 = 9 - Nfam1;
	Nfam2 = 9 - Nfam2;
	Nfam1RI = 9 - Nfam1RI;

	rename Y = RP_Y_SUM
		Nfam0 = RP_Nfam0_SUM
		Nfam1 = RP_Nfam1_SUM
		Nfam1RI = RP_Nfam1RI_SUM
		Nfam2 = RP_Nfam2_SUM;
	RP_N_SUM = sum(of Nfam0 Nfam1 Nfam1RI Nfam2);
	RP_familiar_SUM = sum(of Nfam1 Nfam2); *Nfam1RI;
	RP_SUM = sum(of Y RP_N_SUM);
run;


proc sort data = sumAccurate;
	by subject;
run;

proc sort data = propAccurate;
	by subject;
run;

proc sort data = probesTrans;
	by subject;
run;

data mergedProbes;
	merge probesTrans propAccurate sumAccurate;
	by subject;
run;

proc print data = mergedProbes;
run;

data diss.RecentProbesScores;
	set mergedProbes;
	RP_Familiar_ACC = mean(of RP_Nfam2_ACC RP_Nfam1_ACC);
	if subject in (40, 65, 100) then delete; *near 50% in terms of accuracy;
run;

proc reg data = diss.RecentProbesScores;
	model RP_Familiar_RT = RP_Nfam0_RT;
	output out=probesresid residual=RP_Familiar_RTResid;
quit;

proc sort data = diss.RecentProbesScores;
	by subject;
run;

proc sort data = probesresid;
	by subject;
run;

data diss.RecentProbesScores;
	merge diss.RecentProbesScores probesresid(keep=subject RP_Familiar_RTResid);
	by subject;
run;

proc contents data = diss.RecentProbesScores;
run;

data _null_;
	set diss.RecentProbesScores;
	file "&projectPath\Recent Probes\RecProbes.dat" lrecl=1032;
	put Subject RP_RT RP_Y_RT RP_N_RT RP_Nfam0_RT RP_Nfam1RI_RT RP_Nfam1_RT RP_Nfam2_RT
		RP_ACC RP_Y_ACC RP_N_ACC RP_Nfam0_ACC RP_Nfam1RI_ACC RP_Nfam1_ACC RP_Nfam2_ACC
		RP_SUM RP_Y_SUM RP_N_SUM RP_Nfam0_SUM RP_Nfam1RI_SUM RP_Nfam1_SUM RP_Nfam2_SUM;
run;
