%let projectpath = C:\Documents and Settings\Michael Hallquist\My Documents\Online_Documents\Lab\Dissertation\Analyses;
%let datapath = &projectPath\Data;
*%let datapath = E:\Online_Documents\Lab\Dissertation\Analyses\Data;
libname diss "&datapath";

*setup flanker data for processing;
*note that flanker reaction times max out at 970 (1000ms minus 30ms prerelease);
data flanker;
	set diss.flanker;
	if subject in (7, 14, 36, 38, 47, 55, 62, 83) then delete;
run;

*BEGIN ACCURACY ANALYSES;
proc sort data = flanker;
	by subject;
run;

*following code breaks up parcels by block. Max values stay the same as above;
data altparcels;
	set flanker;
	by subject;
	
	retain C_I_Parcel1 C_I_Parcel2 I_I_Parcel1 I_I_Parcel2;
	retain C_C_Parcel1 C_C_Parcel2 I_C_Parcel1 I_C_Parcel2;

	error = 1-TrialSlide_ACC;

	if first.subject then do;
		C_I_Parcel1 = 0;
		C_I_Parcel2 = 0;
		I_I_Parcel1 = 0;
		I_I_Parcel2 = 0;
		C_C_Parcel1 = 0;
		C_C_Parcel2 = 0;
		I_C_Parcel1 = 0;
		I_C_Parcel2 = 0;
	end;

	if missing(TrialSlide_ACC) then delete;

	*congruent block 1;
	if CongruentTrials_Cycle = 1 then do;
		if Incongruent = 1 then C_I_Parcel1 = C_I_Parcel1 + error;
		else if Incongruent = 0 then C_C_Parcel1 = C_C_Parcel1 + error;
	end;
	else if CongruentTrials_Cycle = 2 then do;
		if Incongruent = 1 then C_I_Parcel2 = C_I_Parcel2 + error;
		else if Incongruent = 0 then C_C_Parcel2 = C_C_Parcel2 + error;
	end;
	else if IncongruentTrials_Cycle = 1 then do;
		if Incongruent = 1 then I_I_Parcel1 = I_I_Parcel1 + error;
		else if Incongruent = 0 then I_C_Parcel1 = I_C_Parcel1 + error;
	end;
	else if IncongruentTrials_Cycle = 2 then do;
		if Incongruent = 1 then I_I_Parcel2 = I_I_Parcel2 + error;
		else if Incongruent = 0 then I_C_Parcel2 = I_C_Parcel2 + error;
	end;

	totalError = sum(of C_I_Parcel1 C_I_Parcel2 I_I_Parcel1 I_I_Parcel2 C_C_Parcel1 C_C_Parcel2 I_C_Parcel1 I_C_Parcel2);
	C_I_Error = sum(of C_I_Parcel1 C_I_Parcel2);
	C_C_Error = sum(of C_C_Parcel1 C_C_Parcel2);
	I_C_Error = sum(of I_C_Parcel1 I_C_Parcel2);
	I_I_Error = sum(of I_I_Parcel1 I_I_Parcel2);
	I_Error = sum(of I_I_Error C_I_Error);
	C_Error = sum(of C_C_Error I_C_Error);
	if last.subject then output;

	keep subject C_I_Parcel1 C_I_Parcel2 I_I_Parcel1 I_I_Parcel2
		C_C_Parcel1 C_C_Parcel2 I_C_Parcel1 I_C_Parcel2
		totalError C_I_Error C_C_Error I_C_Error I_I_Error
		I_Error C_Error;
run;

proc contents data = altparcels;
run;


proc print data = altparcels (where=(subject=1));
	var subject totalError C_I_Parcel1 C_I_Parcel2 I_I_Parcel1 I_I_Parcel2 C_C_Parcel1 C_C_Parcel2 I_C_Parcel1 I_C_Parcel2;
run;

*proc corr data = altparcels;
*	var I_I_Parcel1 I_I_Parcel2 C_I_Parcel1 C_I_Parcel2;
*run;

*proc corr data = altparcels;
*	var I_C_Parcel1 I_C_Parcel2 C_C_Parcel1 C_C_Parcel2;
*run;

proc sql;
	select max(I_I_Parcel1), max(I_I_Parcel2), max(C_I_Parcel1), max(C_I_Parcel2) from altparcels;
quit;

data _null_;
	set altparcels;
	file "&projectPath\Flanker\Parcel.dat" lrecl=1032;
	cerrors = sum(of I_C_Parcel1 I_C_Parcel2 C_C_Parcel1 C_C_Parcel2);
	put I_I_Parcel1 I_I_Parcel2 C_I_Parcel1 C_I_Parcel2 cerrors;
run;

*MERGE RT DATA WITH ACCURACY DATA ABOVE;
proc sort data = altparcels;
	by subject;
run;

data diss.flankerScores;
	merge flankerRTs altparcels;
	by subject;
	C_flankerControl = C_I_MeanRT - C_C_MeanRT;
	I_flankerControl = I_I_MeanRT - I_C_MeanRT;
	flankerControl = I_MeanRT - C_MeanRT;
	if subject = 90 then C_flankerControl = 63.37; *winsorize this score, as it was 121, twice that of anyone else;
run;

/*proc reg data = diss.flankerScores;
	model I_Error = C_Error /r;
	output out=errorScore residual=I_Error_Resid;
quit;*/

proc glimmix data = diss.flankerScores plots=studentpanel(unpack);
	model I_Error = C_Error /solution dist=poisson link=log;
	nloptions technique=nrridg;
	output out=errorScore residual=I_Error_Resid;
run;

proc sort data = errorScore;
	by subject;
run;

data diss.flankerScores;
	merge diss.flankerScores errorScore (keep=subject I_Error_Resid);
	by subject;
run;



proc reg data = diss.flankerscores;
	model I_MeanRT = C_MeanRT;
	output out = residScore residual=I_MeanRT_Resid;
quit;


proc sort data = errorScore;
	by subject;
run;

data diss.flankerScores;
	merge diss.flankerScores residScore (keep=subject I_MeanRT_Resid);
	by subject;
run;

*%inc "&projectPath\CorrWithMacro.sas";
*%corrWith(diss.allcombined, residscore, SNAPFacets);
*%corrWith(diss.allcombined, I_Error, SNAPFacets, partial=C_Error);

data _null_;
	set diss.flankerscores;
	file "&projectPath\Flanker\Parcels.dat" lrecl=1032;
	cerrors = sum(of I_C_Parcel1 I_C_Parcel2 C_C_Parcel1 C_C_Parcel2);
	put Subject 
		I_I_Parcel1 I_I_Parcel2 C_I_Parcel1 C_I_Parcel2 cerrors
		C_I_MeanRT C_C_MeanRT I_C_MeanRT I_I_MeanRT 
		C_I_B1MeanRT C_I_B2MeanRT I_I_B1MeanRT I_I_B2MeanRT 
		C_C_B1MeanRT  C_C_B2MeanRT  I_C_B1MeanRT  I_C_B2MeanRT
	;
run;
