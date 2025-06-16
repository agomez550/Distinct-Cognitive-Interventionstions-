/*Adrian Gomez */ 
/*Final Project Code */
/*Biomethods IV */ 


libname pt 'C:\Users\gomezadr\Documents\Bio4' ;

/*Displaying the data */ 
proc print data=pt.Active (OBS= 6);
run ;

/*Descriptive statistics stratified by intervention group */ 
proc means data=pt.Active N NMISS MIN MAX RANGE mean median std var maxdec=2;
   class INTGRP; 
   var HVLTT HVLTT_2 HVLTT_3 HVLTT_4 HVLTT_5 HVLTT_6;
run;

/*Descriptive statistics for continuous variables */ 
proc means data=pt.Active N NMISS MIN MAX RANGE mean median std maxdec=2; 
   var HVLTT HVLTT_2 HVLTT_3 HVLTT_4 HVLTT_5 HVLTT_6 AGE YRSEDUC MMSETOTL;
run;

/*descriptive statistics for categorical variables */ 
proc freq data=pt.Active;
   TABLES INTGRP DRIVER GENDER SITE; 
run;

/*convtering wide to long form and making baseline its own covariate */ 
data long ; set pt.Active ;
 y=HVLTT;   time=0 ;  baseline= HVLTT; output ;
 y=HVLTT_2; time=1 ;  baseline= HVLTT;  output ;
 y=HVLTT_3; time=12 ; baseline=HVLTT;  output ;
 y=HVLTT_4; time=24 ; baseline= HVLTT; output ;
 y=HVLTT_5; time=36;  baseline= HVLTT; output; 
 y=HVLTT_6; time=60;  baseline=HVLTT; output;
 KEEP t AID HVLTT_2 HVLTT_3 HVLTT_4 HVLTT_5 HVLTT_6 y time MMSETOTL INTGRP REPLCODE SITE BOOSTER GENDER AGE YRSEDUC DRIVER baseline;
run ; 

/*deleted the baseline and treated it as a covariate */ 
data new ; set long; 
if time=0 then delete; 
run; 

/*Code for passing HVLT score (binary response) */ 
data longby2; set new; 
if y ne . then do; 
if y < = 17 then normalscore =0; 
else if y > 17 then normalscore =1; 
end; 
run; 

proc print data=longby2 (OBS= 50);
run ;

/*proportion of passing scores over the 5 time points */ 
proc freq data=longby2;
   TABLES normalscore; 
run;


/*graphing the trend in for individuals */ 
proc sgplot data=long (where=(aid <= 50));
 series x=time y=y / group=AID ;
 yaxis min=0 max=40 label="Mean HVLT Score";
 title " 50 Individual HVLT Intervention Scores "; 
run ;

/* Calculating means by intergroup and time (population) */ 
proc means data=long n mean std nway; 
var y;
class INTGRP time; 
output out= meandata mean=mnHVLTT; 
run; 

/* plotting response profile (population) where baseline is removed */ 
proc gplot data=meandata  (where = (time ne 0)); 
symbol1 color =black
interpol= join
value=dot; 
symbol2 color =blue
interpol =join
value= triangle; 
symbol3 color = red 
interpol = join 
value = circle; 
symbol4 color= green 
interpol= join 
value= square; 
title "Mean Population HVLT Intervention Scores "; 
 legend1 label=("Training Type:")
      value=("Memory Training" "Reasoning Training" "Speed Training" "Control Group")
      position=(bottom center)
      mode=share
	  offset=(0, 5) /* Moves the legend downward */
      frame;

plot mnHVLTT*time=INTGRP / legend=legend1; 
run; 



/*Create quadratic and spline time */ 
data long12; set longby2;    /*set longby2 */ 
timec= time-27; /*(1+12+24+36+60)/5 =26.6=27 */
timesq = time*time ;
timecsq = timec*timec; 
time_1 = max(time-12, 0);
run; 

proc print data=long12 (obs=30); run; 



/*linear mixed model (linear structure)*/
proc mixed data=new method=ml;
 class aid INTGRP  ;
 model y=INTGRP time INTGRP*time / s chisq ;
 random intercept time / type=un subject=aid g s;
run;

/* AIC: 60346.7 -2log-likelihood: 60324.7  */
  

/*quadratic*/
proc mixed data=long12 method =ml;
class AID INTGRP(ref='4'); 
model y=INTGRP timec timecsq timec*INTGRP timecsq*INTGRP/ solution s chisq; 
random intercept timec  / type=un subject= aid g s; 
run; 

/* -2long-likelihood 58742.2 AIC (Smaller is Better) 58774.2 */  
 

/* spline */ 
data spline; 
set new; 
st1= min(time, 12);
st2=max(0, time-12); 
run; 

proc mixed data= spline method =ml; 
class aid INTGRP; 
model y = INTGRP st1 st2 INTGRP*st1 INTGRP*st2 / solution chisq; 
random intercept st1 st2 / type= un subject= aid g s; 
run; 

/* AIC (Smaller is Better) 59428.1, -2 Log Likelihood 59392.1  */ 

/*quadratic with baseline as a covariate */ 



/*quadratic with covariates and baseline */ 
proc mixed data=long12 method =ml;   
class AID INTGRP GENDER(ref='2') SITE DRIVER(ref='0'); 
model y= INTGRP timec timecsq timec*INTGRP timecsq*INTGRP baseline MMSETOTL GENDER AGE YRSEDUC SITE DRIVER /solution s chisq; 
random intercept timec / type=un subject= AID g s; 
run; 
/*AIC: 50040.8, -longlikehood: 49986.8 UPDATE TABLES */ 


/*quadratic with covariates removing baseline */ 
proc mixed data=long12 method= ml;  
class AID INTGRP GENDER(ref='2') SITE DRIVER(ref='0'); 
model y=INTGRP timec timecsq timec*INTGRP timecsq*INTGRP MMSETOTL GENDER AGE YRSEDUC SITE DRIVER / solution s chisq; 
random intercept timec / type= un subject= AID g s;
run; 
/*AIC: 57249.2, -2log likelihood: 57197.2  */ 

proc print data=long12 (obs=30); run; 




/*Passing score status outcomes (including random intercepts*/ 
 
proc glimmix data=long12 method= quad; 
class AID INTGRP(ref="4") GENDER(ref="2") SITE(ref="1"); 
model normalscore =  INTGRP YRSEDUC MMSETOTL GENDER SITE timec INTGRP*timec AGE baseline DRIVER/ d=bin link=logit solution; 
random intercept / subject =AID type=un; 
run; 





