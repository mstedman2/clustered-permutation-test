********************************************************************************;
*   Program:  permutation_macro.sas;
*   Purpose:  To perform a clustered permutation test;
*   input:    /usr/path/dset.sas7bdat;
*   output:   /usr/path/permuteout.sas7bdat;
*   definitions:  datain = input datset;
*                 noperm = number of permutations;
*                 dups = allowed yes /no;
*                 treat = variable for assigned treatment group;
*                 clusterid = cluster id;
*                 seed = seed value for randomly permuting datset;
*                 procout = output dataset from procedure;
*                 testname = variable name of test in procout;
*                 permuteout=name of output dataset for permutation test;
*
*  Copyright (C) 2008  Margaret Stedman;
*  
*    This program is free software: you can redistribute it and/or modify;
*    it under the terms of the GNU General Public License as published by;
*    the Free Software Foundation, either version 3 of the License, or;
*    any later version.;
*
*    This program is distributed in the hope that it will be useful,;
*    but WITHOUT ANY WARRANTY;* without even the implied warranty of;
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the;
*    GNU General Public License for more details:;
*    <http://www.gnu.org/licenses/>;
******************************************************************************;
options nofmterr;

******************************************************************************; 
******Enter directory for input dataset here:*********************************;

libname ll "C:\My SAS Files\9.1\";

******Enter directory for output dataset here:********************************;

libname ll2 "C:\My SAS Files\9.1\";

******Enter procedure within macro statements here:***************************;

%macro procmacro;
%mend procmacro;

******************************************************************************;
******************************************************************************;

 %macro ptest (datain, noperm, dups, treat, clusterid, seed, procout, testname,
  permuteout);

data datain;
  set ll.&datain;
run;

title "Original Result";
%procmacro;

proc sql;
   title "Number of rows in output dataset, &procout";
   select count(*) into: lines
   from &procout;

%if %eval(&lines) ne 1 %then %do;
   %put ERROR: output dataset, &procout, has more than 1 row;
%end;

   create table ll.&permuteout as 
   select *, 0 as PLoup 
   from &procout; 

   title "Number of clusters in dataset";
   select count (distinct &clusterid) into : noclus 
   from datain;

   title "Number of clusters per group";
   select count (distinct &clusterid) into :noclus1 - :noclus2
   from datain
   group by &treat;

   title "Number of possible permutations";
   select comb(&noclus, &noclus1) into :block
   from &procout;

quit;

%if &block <= &noperm %then &noperm = all;

proc sort data=datain;
   by &clusterid;
run;

data dat; 
    set datain; 
    by &clusterid; 
    if _N_=1 then clusterid=0;
    retain clusterid; 
    if first.&clusterid then clusterid=clusterid +1; 
	if last.&clusterid then output;
     keep &clusterid clusterid;
run; 

%if &noperm ^= all %then %do;
proc iml; 
x=j(&noclus, &noperm,.); 
g=j(&noclus, &noperm,.); 
call randseed (&seed); 
call randgen(x,'uniform');
%do i=1 %to &noperm;
call sortndx(ndx,x,{&i});  
g[,&i]=choose(ndx[,1] > &noclus1 ,1,0); 
%end; 
gt=g`;
create outdat from gt; 
append from gt; 
quit; 

%if &dups = no %then %do;
proc sort data	= outdat nodup;
 by _all_;
 run;
%end;

proc sql;
   title "Number of permutations performed"; 
   select count(*) into : noperm2
   from outdat;
quit;

 proc transpose data=outdat out=outdat2;
 run;
 
%end;
 

%else %if &noperm = all %then %do;

data fullset;
   do col=1 to &block;
   do clusterid=1 to &noclus;
   output;
   end;
   end;
run;

proc plan seed=&seed;
factors col = &block
ordered
clusterid = &noclus1 of &noclus comb /
noprint;
output out = combx;
run; quit;

data idcombx;
  merge fullset (in=a) combx (in=b);
  by col clusterid;
  treat=0;
  if b then treat=1;
run;

proc transpose data=idcombx out=combxt;
   var treat;
   by col;
run;
quit;

proc transpose data=combxt (drop= col _name_) out=outdat2;
run;

%let noperm2 = &block;

%end;

proc sql; 
  create table idcluster as
  select a.*, b.&clusterid
  from outdat2 a, dat b
  where input(compress(a._name_,"COL"),8.) = b.clusterid; 
quit; 

data notreat;
   set datain;
   drop &treat;
run;

%do p=1 %to &noperm2; 

proc sql;
create table analdat as 
   select b.col&p as &treat, a.* 
   from notreat a, idcluster b 
   where a.&clusterid = b.&clusterid; 
quit;
  
ods listing close;
 
%procmacro;

data ll.&permuteout; 
    set ll.&permuteout (in=a) &procout (in=b); 
    if b then PLoup=&p; 
run; 
 
%end; 

proc sql;
   create table origout as
   select *
   from ll.&permuteout
   where PLoup=0;

   create table permuted as
   select * 
   from ll.&permuteout
   where PLoup>0;
 
   create table permuteout as 
   select a.*, b.&testname as obstest, 
   case 
     when abs(obstest)le abs(a.&testname) then 1
     else 0 
   end as test format=8.
   from permuted a, origout b;
quit;

ods listing;
proc means data=permuteout noprint;
    var test;
	output out=results mean=p_value n=permutations;
 run;
 
 proc print data=results noobs;
 var permutations p_value;
 Title "Significance testing for Permutation test";
 format p_value 8.7;
 run;

%mend ptest;
quit;

***Macro call;
%ptest (datain=, noperm=, dups=, treat=, clusterid=, seed=, procout=, testname=,
 permuteout=);

 /*
***Example 1:  Chi-square test;
proc freq;
  table bmd_op_fol * groupn/chisq;
  output out=procout chisq;
run;

%ptest (datain=horizdat, noperm=1000 , dups=yes, treat=groupn, clusterid=dr_id,
 seed=34637, procout=procout, testname=_PCHI_, permuteout=result1);
 
***Example 2:  Logistic Regression;
ods output parameterestimates=procout (where=(variable="groupn")); 
proc logistic descending;
   model bmd_op_fol = groupn;
run;

%ptest (datain=horizdat, noperm=1000 , dups=yes, treat=groupn, clusterid=dr_id,
 seed=345, procout=procout, testname=WaldChiSq, permuteout=result1);
*/
