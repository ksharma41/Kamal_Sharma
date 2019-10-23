helix30.c
Long ago
Mar 3, 2017
K
You uploaded an item
C
helix30.c
#include <stdio.h>                    /*think about relation between random numbers and modify program if necessary*/
#include <stdlib.h> 
#include <time.h>
#include <math.h>


#define STEP 0.261               /*starting step length in radians*/
#define FIELDSTART 0.0                                
#define FIELDSTOP 2.533333e-6
#define FIELDSTEP 1.26e-7


#define CELCIUS 173.5                /*Actual Physical Temperature of System*/ 
#define TAF 175.0
#define TPU 187.0

#define A 1.3
#define B 13.0
#define J2 1.0e-5
#define j 5.0e-5
#define LAYERS 60					
#define TINIT 10.0              /*Initial dummy temp*/
#define TFINAL 1.0e-17              /*Final dummy temp*/
#define TFACTOR 0.875              /*Temperature Reduction Coefficient*/

#define TLOOP 300       /*Maximum number of cycles to be carried out at a particular temperature*/
#define NS 20                /*no of cycles for a particular step length*/
#define D 2.0                   /*factor by which step to be multiplied if success rate is 100% in last cycle*/
int GetRand(int min, int max);
int GetRandblah(int min, int max);


main()
{

FILE *gp,*ep,*xp;
gp = fopen("Complete Configuration.txt","w");
ep = fopen("Energy_Field.txt","w");
xp = fopen("Field vs Polarization.txt","w");
time_t start,end;


/*variable declaration*/     
  int u,g,i,q,r,b, layer=LAYERS, tindex, k,  success[LAYERS],moves[LAYERS];
  double runtime,fraction,vnew[layer],vold[layer];
  double J1,PI,J3,c,temp=TINIT,initial[layer],step[layer],evolve_from[layer],evolve_to[layer],optimum[layer],globalmin[layer];
  double field,x,a,PSQRD,ediff,boltzman,initvalue,newvalue,oldvalue,optvalue,globmin=9999.9,ppl;
start=time(NULL);
PI=acos(-1);
J1=j*(CELCIUS-TAF);
PSQRD=A*(TPU-CELCIUS)/B;
field=FIELDSTART;
J3=(2*J2*(PSQRD,2)-J1*PSQRD/cos(30*PI/180))/4;


/*  fprintf(fp,"*********************NEW RUN*********************  \n \n \n \n");
  */ 

fprintf(xp,"Numbers of Layers=%d\nField\t\tPolarization Per Layer\n",(LAYERS-10));
                
decfield:
         temp=TINIT;
         initvalue=newvalue=oldvalue=optvalue=0; globmin=9999.9;
         
             for (i = 0; i < layer; i++)
  { 
    r = GetRand(0, 1000);                       /*change this 1000 to any comfortable value depending on accuracy*/
    initial[i]=(double) r*PI*2/1000;            /*initial random state established between 0 to 2pi*/   
    //fprintf(fp,"Initial angle %d is %f \n",i,initial[i]*180/PI);
    vnew[i]=STEP;
    vold[i]=0;
    moves[i]=1;}          
   
for (i = 0,c=0,initvalue=0; i < layer-2; i++)
{
    c=J1*PSQRD*cos(initial[i+1]-initial[i]) - J2*pow(PSQRD,2)*pow(cos(initial[i+1]-initial[i]),2) + J3*cos(initial[i+2]-initial[i]) - sqrt(PSQRD)*field*cos(initial[i]);                                          
initvalue=initvalue+c;     /*calc of variable part of free energy for initial random config*/
        }
oldvalue=initvalue;    /*store initial energy into oldvalue*/


/*Assigning the initial values to a new array in the i dimensional space*/
for(i = 0; i < layer; i++)
{
      evolve_from[i]=initial[i];                     /*leaving initial point undisturbed for experimental purposes*/
      /*fprintf(fp,"%f \n",evolve_from[i]*180/PI);
      */
      }

/********************************************/
/*********START OF ANNEALING*****************/
/********************************************/

anneal:
       
       if(temp<TFINAL) goto finish; 
       
       /*fprintf(fp,"Temperature of the system : %f \n",temp);
       */
       tindex=1;
temprepeat:

if(tindex<=TLOOP )  /*checks modulus of change in energy with ETOL which i removed lalter on, its just temp check now*/
{                 

/*printf("Cycle Number # %d \n",tindex);
fprintf(fp,"Cycle Number # %d \n",tindex);*/

           /*generating random variable to get random displacement from evolve_from position*/

for(q=0;q<layer;q++) {moves[q]=1; success[q]=1;}


for(q=0;q<=NS;q++)
{
u=0;
nextdimension:
//fprintf(fp,"***********************\nstep = %f \n",vnew[u]*180/PI);
if(u<layer)
{  /*this paranthesis end before tindex thingy*/
{
    r = GetRand(0,1000);                     /*change this 2000 to any comfortable value depending on accuracy*/
    fraction = (double) -1 + (double) r/500;              
    evolve_to[u]=evolve_from[u]+vnew[u]*fraction;
    evolve_to[u]=evolve_to[u]-(floor(evolve_to[u]/(2*PI))*2*PI);
   
 } 

/*Computing value of the energy function in evolve_to configuration*/
for(i = 0,newvalue=0; i < layer-2; i++)
{    
      a=J1*PSQRD*cos(evolve_to[i+1]-evolve_to[i]) - J2*pow(PSQRD,2)*pow(cos(evolve_to[i+1]-evolve_to[i]),2) + J3*cos(evolve_to[i+2]-evolve_to[i]) - sqrt(PSQRD)*field*cos(evolve_to[i]-PI);
      newvalue=newvalue+a;   /*calculating new value of free energy for proposed state*/
      }
//newvalue=newvalue+ J1*PSQRD*cos(evolve_to[layer-1]-evolve_to[layer-2]) - J2*pow(PSQRD,2)*pow(cos(evolve_to[layer-1]-evolve_to[layer-2]),2) + J3*cos(evolve_to[0]-evolve_to[layer-2]) - sqrt(PSQRD)*field*cos(evolve_to[layer-2]);     
//newvalue=newvalue+ J1*PSQRD*cos(evolve_to[0]-evolve_to[layer-1]) - J2*pow(PSQRD,2)*pow(cos(evolve_to[0]-evolve_to[layer-1]),2) + J3*cos(evolve_to[1]-evolve_to[layer-1]) - sqrt(PSQRD)*field*cos(evolve_to[layer-1]);
   
if(newvalue<oldvalue)
{
/*fprintf(fp,"New Config Accepted");*/
success[u]=success[u]+1;
moves[u]=moves[u]+1;/*fprintf(fp,"successes for layer no %d %d\n",u,success[u]);
*/optvalue=newvalue;           /*storing the energy in optvalue*/
for(i = 0; i < layer; i++)
{
      optimum[i]=evolve_to[i];     /*storing the stable config in optimum[]*/  
} /* fprintf(fp,"$$ New Angle %d is %f\t",u,evolve_to[u]*180/PI);
printf("Free Energy of above configuration is %E \t\t ***** \n",optvalue);
fprintf(fp,"lower energy = %f\n",optvalue); 
*/}
else
{
    ediff=oldvalue-newvalue;
    boltzman=exp(ediff/temp);      /*probability of system to be found in proposed configuration*/
   // fprintf(fp,"Boltzman = %f\n",boltzman);
    if(boltzman>=0.9) boltzman=0.9;
    b=GetRand(0,1000);
    x=(double) b/1000;  
moves[u]=moves[u]+1;
    
    if(x<boltzman)
{  //fprintf(fp,"Boltzman = %f\n",boltzman);
    
    success[u]=success[u]+1;
    optvalue=newvalue;
    /*fprintf(fp,"@@ New Angle %d is %f \t",u,evolve_to[u]*180/PI);*/
    for(i = 0; i < layer; i++)
      optimum[i]=evolve_to[i];  
//printf("Free Energy of above configuration is %E \n",optvalue);
/*fprintf(fp,"HIGHER ENERGY = %f\n",optvalue);*/
}
else
{
//printf("New Random Configuration Rejected on acount of Boltzmann Probability Distribution \n");
/*fprintf(fp,"New Random Configuration Rejected on acount of Boltzmann Probability Distribution \n");
*/
}

}
for(i = 0; i < layer; i++)
{
  //    printf("temp =%e, loop = %d, layer = %d\n %d\n",temp,q,u,tindex);
      evolve_from[i]=optimum[i];     /*storing the stable config in optimum[]*/  
} 
oldvalue=optvalue;

//fprintf(fp,"Cycle Number %d For Layer %d\tTotal = %d\tAccepted=%d  with Probability =%f\n",tindex,u,moves[u],success[u],success[u]/moves[u]);
u=u+1;
goto nextdimension; 
}
} /*end of biggest FOR with q varuiable*/

/*STEP VARIATION-entire step vector at a time*/
for(g=0;g<layer;g++)
{
vold[g]=vnew[g];
if(success[g]>0.6*NS)  vnew[g]=vold[g]*(1+D*((success[g]/NS)-0.6)/0.4);
if(success[g]<0.4*NS)  vnew[g]=vold[g]/(1+D*(0.4-(success[g]/NS))/0.4); 
vnew[g]=vnew[g]-(floor(vnew[g]/(2*PI))*2*PI); /*gives remainder of new step when divided by PI */
 
if(vnew[g]>PI) vnew[g]=2*PI-vnew[g]; 
/*fprintf(fp," Step for Angle %d is %f\n",g,vnew[g]*180/PI);
 /*fprintf(fp,"Success %d\n",success[g]);
//fprintf(fp,"%f\n",vnew[g]*180/PI);*/
}

tindex=tindex+1;

goto temprepeat;
}


else { 
     
     temp=TFACTOR*temp; 

end=time(NULL);
runtime=difftime(end,start);
 


goto anneal;}
/*printing the final answer*/
finish:

if(optvalue<globmin)
{
for(i = 0; i < layer; i++)
globalmin[i]=optimum[i];
globmin=optvalue;
}
 
fprintf(gp,"\n***********************\n***********************\n");
fprintf(gp,"Next Stable Configuration at FIELD=%e  Physical Temperature = %f\n is as follows:-\n",field,CELCIUS); 


for(i = 0; i < layer; i++)
{fprintf(gp,"%f\n",(globalmin[i]*180/PI)); 
}
 
fprintf(gp,"Free Energy=%e Joules\t",globmin);
fprintf(gp,"Electric Field=%e \n \n \n \n",field);
end=time(NULL);
runtime=difftime(end,start);
fprintf(ep,"%e\t%e\n",globmin,field);

for(i = 5,ppl=0; i < layer-5; i++)
{
ppl=ppl+sqrt(PSQRD)*cos(globalmin[i]);
}
fprintf(xp,"%e\t%e\n ",field,ppl/(layer-10)); 


field=field+FIELDSTEP;
if(field<=FIELDSTOP)


{
goto decfield;
}

fprintf(gp,"Run Duration: %f seconds\nJ1=%f  J2=%f  J3=%f P=%f *****************************\n**************************\n END OF RUN\n",runtime,J1,J2,J3,sqrt(PSQRD)); 
fclose(gp); 
fclose(ep); 
fclose(xp);
return(0);
}






int GetRand(int min, int max)
{ 
  static int Init1 = 0;
  int rc;
  
  if (Init1 == 0)
  {
    /*
     *  As Init is static, it will remember it's value between
     *  function calls.  We only want srand() run once, so this
     *  is a simple way to ensure that happens.
     */
    srand(time(NULL));
    Init1 = 1;
  }

  /*
   * Formula:  
   *    rand() % N   <- To get a number between 0 - N-1
   *    Then add the result to min, giving you 
   *    a random number between min - max.
   */  
  rc = (rand() % (max - min + 1) + min);
  
  return (rc);
}

int GetRandblah(int min, int max)
{ 
  static int Init = 0;
  int rc;
  
  if (Init == 0)
  {
    /*
     *  As Init is static, it will remember it's value between
     *  function calls.  We only want srand() run once, so this
     *  is a simple way to ensure that happens.
     */
    srand(time(NULL));
    Init = 1;
  }

  /*
   * Formula:  
   *    rand() % N   <- To get a number between 0 - N-1
   *    Then add the result to min, giving you 
   *    a random number between min - max.
   */  
  rc = (rand() % (max - min + 1) + min);
  
  return (rc);
}