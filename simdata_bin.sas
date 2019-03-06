********************************************************************************;
*   Program:  sim_data.sas;
*   Purpose:  To simulate data for example study;
*   output:   /usr/path/simdat.sas7bdat;
*   definitions:  b00=  intercept parameter;
*				  b01=  treatment parameter;
*				  tau=  variance of random effect;
*				nclust= number of clusters;
*				ntreat= number of treated clusters;
*				nsub= number of subjects per cluster;
*				seed= seed for simulation;
******************************************************************************;
options nofmterr;

******************************************************************************; 
******Enter directory for input dataset here:*********************************;

libname ll "N:\Projects\Other or mixed\Margaret\clustered logrank";

%macro sim(b00,b01,tau,nclust,ntreat,nsub,seed);
data simdatb;
   informat id $7.;
      sseed=&seed;
      do j=1 to &nclust;
      if j<=&ntreat then treat=1; /****treatment group********/
      else treat=0;
      clustid=j; /*******cluster id***********/
      call rannor(sseed,z);
      u0j=z*&tau; /***********random intercept************/
        do i=1 to &nsub;
           idp1=put(j,$3.);
           idp2=put(i,$3.);
           id=compress(idp1||'-'||idp2); /********subject id*********/
           mod=&b00+&b01*treat+u0j;
           p=exp(mod)/(1+exp(mod));
           call ranbin(sseed,1,p,event);
          output;
        end;
    end;
run;

%mend sim;
quit;

%sim(b00=-2.3,b01=.37,tau=1,nclust=458,ntreat=229,nsub=8,seed=656);

data ll.simdatb;
   set simdatb;
run;

proc contents data=ll.simdatb;
run;
