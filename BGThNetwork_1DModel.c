/****************************************************************************************************************

****************************************************************************************************************/
/**
* @file BGThNetwork_1DModel.c
* @author Osvaldo M Velarde
* @date May 2017
* @brief Simulation of a Cortico-Basal Ganglia-Thalamocortical network model.
* The architecture and dynamics of the network can be seen in the reference.
* The objective is to analyze the effect of the intensities of the synaptic connections on the dynamics of the network.
* In particular, associating pathological and physiological states in Parkinson's disease.
* Furthermore, it is desired to study the effect of deep brain stimulation on the pathological dynamics of the network.
* @see https://doi.org/10.1371/journal.pone.0182884
*/

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include<time.h>

/* Simulation times (ms)*/
#define t_0_sim 0.0	        // Start time
#define t_f_sim 6000.0      // End time
#define	dt 0.5 				// Temporal step 

/* External input in Cortex (Ctx). Time interval = [t_0_Ctx, t_f_Ctx] (ms)*/
#define t_0_Ctx 250.0
#define t_f_Ctx 6000.0

/* External input in Striatum (Str). Time intervals = [t_0_Str1, t_f_Str1] U [t_0_Str2, t_f_Str2](ms)*/
#define t_0_Str1 300.0
#define t_f_Str1 1000.0
#define t_0_Str2 3000.0
#define t_f_Str2 4000.0

/* Deep Brain Stimulation (DBS). Time interval = [t_0_DBS, t_f_DBS] (ms)*/
#define t_0_DBS 1000.0
#define t_f_DBS 6000.0

#define PI  3.1415926

/*------------------ Structures -------------------------*/
/*-------------------------------------------------------*/

/**
* @brief neuron_struct represents to a neuron in a population.
*/

typedef struct neuron_struct{
    double parameters[3];                               /**< Values of "threshold", "external input", "noise level". */
    double input_det;                                   /**< Total input. */
    double input_sto;                                   /**< Noise. */
    double angle;                                       /**< Position. */
    double (*ext)(double ext_p[5],double k,double t);   /**< External input. */
    }neuron;

/**
* @brief BGan_struct represents to a neuronal population.
*/

typedef struct BGan_struct{
    neuron *circle;            /**< Array of neurons. */
    int neuron_number;         /**< Number of neurons in a circle. */
    }BGan;                       

/**
* @brief connection_struct represents to a synapse between two neurons (pre and post- synaptic).
*/

typedef struct connection_struct{
    double sinaptic_curr;           /**< Synaptic current m_ij. */
    struct neuron_struct *last;     /**< Pre-synaptic neuron i. */    
    struct neuron_struct *next;     /**< Post-synaptic neuron j. */
    }connection;                 

/**
* @brief cable represents to an array of synapses. All synapses in the "cable" have the same post-synaptic neuron.
*/
typedef connection* cable;

/**
* @brief conductor saves the information on the parameters of the connections between populations.
*/

typedef struct {
    double parameters[3];           /**< Values of synaptic efficacy, time constant and delay. */
    cable* cond;                    /**< Array of "cable". Each element of the array is asociated with a neuron in post-sinaptic population*/
    BGan* last;                     /**< Pre-synaptic Population.  */
    BGan* next;                     /**< Post-synaptic Population. */ 
	int K;							/**< Average number of connections per neuron. */
	double level_div;				/**< Variance of number of connections. */
	}conductor;


/*------------------ Functions: Use of angles --------------------------------*/
/*----------------------------------------------------------------------------*/

/**
* @brief Wrap a double to the interval [-pi,pi) 
* @param value Angle (rad).
* @return Angle (rad) between [-pi, pi].
*/

double normalize_angle(double value){               
    double v1=cos(value);
    double v2=sin(value);
    double theta=acos(v1);
    if(v2>=0)
        return theta;
    else
        return -theta;}

/**
* @brief Sum of angles 
* @param theta1 Angle (rad) in [-pi,pi).
* @param theta2 Angle (rad) in [-pi,pi).
* @return Sum of theta1 and theta2 (rad) in [-pi,pi).
*/

double sum_angle(double theta1,double theta2){
    theta1=normalize_angle(theta1);
    theta2=normalize_angle(theta2);
    double theta= theta1+theta2;
    theta=normalize_angle(theta);
    return theta;}


/*----------------- Connectivity: Probability of connection ------------------*/
/*----------------------------------------------------------------------------*/
/**
* @brief Function required to calculate the probability of connection between two neurons.
* @param x Variable
* @param sigma Parameter of function.
* @return exp((x-1.0)/sigma^2)
*/

double function_gauss(double x,double sigma){
    double e=exp((x-1.0)/(sigma*sigma));
    return e;
}

/**
* @brief Approximation of Gaussian Function in angle variables. 
* The probability of connection between two neurons depends of the angular distance.
* d_12 = 1 - Cos(theta1 - theta2)
* Probability is proporcional to exp(d_12/sigma^2)
* @param theta1 Angle (rad) in [-pi,pi)
* @param theta2 Angle (rad) in [-pi,pi)
* @param sigma Variance.
* @return Probability of connection between neuron 1 (in theta1) and neuron 2 (in theta2)
*/

double probability_connect(double theta_1,double theta_2,double sigma){
    double delta=sum_angle(theta_1,-theta_2);
    double x=cos(delta);
    double e=function_gauss(x,sigma);
    return e;
}

/**
* @brief Computer of normalization constant for the distribution of connections-
* @param A Pre-synaptic population 
* @param i Index for neuron in A. Integer value in [0, Num_neurons_A)
* @param B Post-synaptic population 
* @param K_ab Average number of connections per neuron "i". 
* @param sigma Variance.
* @return Normalization constant
*/

double cte_normalizat(BGan A,int i,BGan B,int K_ab,double sigma){
    int j;
    double No=0;
    K_ab=1.0*K_ab;
    for(j=0;j<B.neuron_number;j++){
         No=No+ probability_connect(A.circle[i].angle,B.circle[j].angle,sigma);
        }
    No=No/K_ab;
    return No;
}

/*------------------ Electrode: Spatial Modulation --------------------------*/
/*---------------------------------------------------------------------------*/
/**
* @brief Spatial modulation for measurement and stimulation by electrodes.
* @param position Position of electrode [-pi,pi) 
* @param threshold Threshold of electrode [0,1)
* @param vol_tissue Angular volume of tissue involved in measurement or stimulation (rad) 
* @param angle Position of neuron.  
* @return Spatial modulation. 
*/

double spatial_electrode(double position, double threshold, double vol_tissue, double angle){
    double z=(cos(angle-position)-1)/(cos(vol_tissue/2)-1);
    z=(1.0-1.0/threshold)*z;
    return 1.0/(1-z);
}

/*----------------------------- NOISE & INPUTS--------------------------------*/
/*----------------------------------------------------------------------------*/
/**
* @brief Generation of random numbers (normal distribution)
* @param mean Mean of distribution position Position of electrode [-pi,pi) 
* @param sigma Standard desviation of distribution
* @return Random number (White noise).
*/

double external_noise(double mean,double sigma){
    double x,y,z;
    x=(double)rand()/(double)RAND_MAX;
    y=(double)rand()/(double)RAND_MAX;

    if(x!=0 && y!=0){
       z=sqrt(-2.0*log(x))*cos(4.0*y*asin(1.0));
       z=sigma*z+mean;
       return z;}


    if (x==0 && y!=0){
       z=sqrt(-2.0*log(y));
       z=sigma*z+mean;
       return z;}

     if (x!=0 && y==0){
       z=sqrt(-2.0*log(x));
       z=sigma*z+mean;
       return z;}

    else
        return 0;
}

/** OBSERVATION: All external_inputs functions follow the form:
* (a + b*cos(angle)) * (c + d*cos(2pi *t/T) )
* def_parameters = [a,b,c,d,T]
* angle: position of neuron
* t: time.
*/ 


double external_default(double def_parameters[5],double angle, double t){
		return def_parameters[0];
}

double external_Ctx(double Ctx_parameters[5],double angle, double t){
        if (t_0_Ctx<t && t<t_f_Ctx){
        double spatial_comp=Ctx_parameters[0]+Ctx_parameters[1]*cos(angle);
        double temporal_comp=Ctx_parameters[2]+Ctx_parameters[3]*cos(2*PI*t/Ctx_parameters[4]);
        return spatial_comp*temporal_comp;
        }
		return 0;
}

double external_Str(double Str_parameters[5],double angle, double t){
        if ((t_0_Str1<t && t<t_f_Str1)||(t_0_Str2<t && t<t_f_Str2)){
        double spatial_comp=Str_parameters[0]+Str_parameters[1]*cos(angle);
        double temporal_comp=Str_parameters[2]+Str_parameters[3]*cos(2*PI*t/Str_parameters[4]);
        return spatial_comp*temporal_comp;
        }
		return 0;
}

/*---------------------- Deep Brain Stimulation ---------------------------*/
/*-------------------------------------------------------------------------*/
double DBS_Period(double t,double amplitude,double width, double slope){
    if (0<=t && t<slope)
        return amplitude*t/slope;
    if(slope<=t && t<=slope+width)
        return amplitude;
    if(slope+width<t && t<2*slope+width)
        return amplitude*(2 - (t-width)/slope);
    else
        return 0;
}

double DBS_Periodic(double t, double amplitude, double period,double width,double slope){
    if (t_0_DBS<=t && t<t_f_DBS){
            double theta=2*PI*(t-t_0_DBS)/period;
            double y=acos(cos(theta));
            if (sin(theta)>0)
                return DBS_Period(period*(y/(2*PI)),amplitude,width,slope);
            else
                return DBS_Period(period*(1-y/(2*PI)),amplitude,width,slope);}
    else
        return 0.0;
}

double random_frequency(double mean_frequency,double coef_var){
    if (coef_var==0)
        return mean_frequency;
    else{
        double p=1/pow(coef_var,2);
        double a=p/mean_frequency;
        double x_1=rand()/(double)RAND_MAX;
        double x_2=rand()/(double)RAND_MAX;
        double e_1=-log(1-x_1);
        double e_2=-log(1-x_2);
        while(e_2 < (p-1)*(e_1- 1 -log(e_1))){
            x_1=rand()/(double)RAND_MAX;
            x_2=rand()/(double)RAND_MAX;
            e_1=-log(1-x_1);
            e_2=-log(1-x_2);}
        double z=p*e_1/a;
        return z;}
}

double DBS_Random(double t,double amplitude, double width,double slope, double* start,int number_pulse){
    int i;
    for(i=0;i<number_pulse;i++){
        if (start[i]<=t && t<=start[i]+width)
            return amplitude;
        if (start[i]+width<t && t<start[i+1])
            return 0;
    }
    return 0;
}

/*------------------- Network construction -----------------------------------*/
/*----------------------------------------------------------------------------*/
int Network_construction(BGan **Population,int Number_BG,conductor **Network,int ***Numbers_Connection){

	int i,j,k; 						// Auxiliar indexs
	double p1,p2,p3; 				// Parameters
	int N_neuron, N_conductor;		// Number of neurons 
	int a,b; 						// Last,Next BGanglia (Conductor)
	int K,sigma;					// Nº average of connections and level of divergence 
	int n,m;						// Number of neurons in last and next BG
	int K_max;						// Max number of connections (aproximmate)
	int count;						// Number of connections in each cable.
	double No;						// Constant of normalization
	double p,x;						// Probability.

	//Data files: Conductors, BGanglia and propierties	
	FILE *Interaction; 
	FILE *Parameters_C;
	FILE *Parameters_G;

	//Create Number_BG
	*Population= (BGan*)malloc(Number_BG*sizeof(BGan));

	//Create neurons in each BGan - Load parameters.
	Parameters_G=fopen("Ganglios.dat","r");
	for(i=0;i<Number_BG;i++){
    	fscanf(Parameters_G,"%lf \t %lf \t %lf \t %d\n",&p1,&p2,&p3,&N_neuron); // Threshold - External input - Noise - Number_neurons
    	(*Population)[i].neuron_number=N_neuron;
    	(*Population)[i].circle=(neuron*)malloc(N_neuron*sizeof(neuron));
        	for(j=0;j<N_neuron;j++){
        		(*Population)[i].circle[j].angle=-PI+(j*2*PI)/(1.0*N_neuron);
        		(*Population)[i].circle[j].parameters[0]=p1;
        		(*Population)[i].circle[j].parameters[1]=p2;
        		(*Population)[i].circle[j].parameters[2]=p3;
        		(*Population)[i].circle[j].input_det=0;
        		(*Population)[i].circle[j].input_sto=0;
        		(*Population)[i].circle[j].ext=external_default;
        	}
	}
	
	fclose(Parameters_G);

	//Conductors between BGanglia
	Interaction=fopen("Interaction.dat","r");
	i=0;
	
	while (feof(Interaction)==0){
		fscanf(Interaction,"%d\n",&a);
		i=i+1;
	}

	N_conductor=i/4;
	rewind(Interaction);

	//Create conductors.
	*Network= (conductor*)malloc(N_conductor*sizeof(conductor));

	//Create a matrix for count number of connections.
	*Numbers_Connection=(int**)malloc(N_conductor*sizeof(int*));

	//For each conductor, create connections
	for(i=0;i<N_conductor;i++){
    
		//Last BGanglia and Next BGanglia
    	fscanf(Interaction,"%d \t %d \t %d \t %d \n",&a,&b,&K,&sigma);
    	(*Network)[i].last=&(*Population)[a-1];
    	(*Network)[i].next=&(*Population)[b-1];
    	n=(*Network)[i].last->neuron_number;
    	m=(*Network)[i].next->neuron_number;
    
		//Struct propierties
		(*Network)[i].K=K;
		(*Network)[i].level_div=1.0*sigma/100;
	
	    //For each conductor, create "n" (last_neuron_number) cables
    	(*Network)[i].cond=(cable*)malloc(n*sizeof(cable));
    	(*Numbers_Connection)[i]=(int*)malloc(n*sizeof(int));

	    //Maximum number of possible connections (aprox).
    	K_max=K+10*floor(sqrt(1.0*K));

	    for(j=0;j<n;j++){
        
			//In each cable, there are K connections (aprox). Upper bounded is K_max.
        	((*Network)[i].cond)[j]=(connection*)malloc(K_max*sizeof(connection));
			count=0;    
        	No=cte_normalizat(*((*Network)[i].last),j,*((*Network)[i].next),(*Network)[i].K,(*Network)[i].level_div);  //NORMALIZACION.

        	//Using random laws, create connections.
        	for(k=0;k<m;k++){
            	x=rand()/(double)RAND_MAX;
            	p=probability_connect((*(*Network)[i].last).circle[j].angle,(*(*Network)[i].next).circle[k].angle,(*Network)[i].level_div)/No;
            	if(x<p){
                	((*Network)[i].cond)[j][count].last=&(*Population)[a-1].circle[j];
                	((*Network)[i].cond)[j][count].next=&(*Population)[b-1].circle[k];
                	count=count+1;}
        	}
    		
    		(*Numbers_Connection)[i][j]=count;
    	}
	}
	
	fclose(Interaction);

	//Load parameters of connections.
	Parameters_C=fopen("Conductores.dat","r");
	for (i=0;i<N_conductor;i++){
    	fscanf(Parameters_C,"%lf \t %lf \t %lf\n",&p1,&p2,&p3);
    	(*Network)[i].parameters[0]=p1;
    	(*Network)[i].parameters[1]=p2;
    	(*Network)[i].parameters[2]=p3;
	}
	fclose(Parameters_C);

	return N_conductor;
}

/*-----------------------LISTAS COMO VECTORES---------------------------------*/
/*----------------------------------------------------------------------------*/
int insertarfinal(double *v, int *pnum, double value,int MAX) {
  if((*pnum)==MAX)
    return -1;

  v[*pnum]=value;
  (*pnum)++;
  return 0;}

int borrarinicio(double *v, int *pnum, double *pvalue) {
  int i;

  if((*pnum)==0)
    return -1;

  *pvalue = v[0];
  for(i=1;i<(*pnum);i++)
    v[i-1]=v[i];

  (*pnum)--;
  return 0;
  }


/*------------------- CAMBIAR INPUT EXTERNO-----------------------------------*/
/*----------------------------------------------------------------------------*/
void Change_InputExt(BGan **Population,double (*New_function)(double ext_param[5],double k, double t),int index_BG, int index_neuron){
    if (0<=index_neuron && index_neuron<(*Population)[index_BG-1].neuron_number){
        (*Population)[index_BG-1].circle[index_neuron].ext=New_function;
    }
    else
        printf("Error al Cambio %d\n",index_neuron);
}

/*------------------- FUNCION DE ACTIVACION- ---------------------------------*/
/*----------------------------------------------------------------------------*/
double Activation_function(double det_input,double stoch_input,double threshold,double *aux){
    double z;
    double total_input=det_input+stoch_input;
    if (total_input<threshold){
        *aux=0.0;
        return 0.0;}
    else
        {z=det_input-threshold;
        *aux=stoch_input;
        return z;}
}

void solve_equations(BGan *Population,int Number_BG,conductor *Network,int N_BGconnection, double ****data_delayed,int dim_delay[],int **Numbers_Connection,int BG_Spike,double **ext_parameters,double **LFP_parameter,double *Electr_parameter,double *DBS_parameter){
//void solve_equations(BGan *Population,int Number_BG,conductor *Network,int N_BGconnection, double ****data_delayed,int dim_delay[],int **Numbers_Connection,int BG_Spike,double **ext_parameters,double **LFP_parameter,double *Electr_parameter,double *DBS_parameter,double Variability){

	//Inputs:
	//-Population[Number_BG]
	//-Network[N_BGconnection]
	//data_delayed[BG_index,neuron_index,connection,time]
	//dim_delay[N_BGconnection]
	//Numbers_connection[BG_index, neuron_index]
	//index BG_spike
	//[H_Ctx,H_str] LFP_parameters(g) DBS_parameters(g) DBS_parameters(t) variability

	int i,j,k,s;   		//indexs
	int N_aux,Kab;      //aux integers

	double value, aux; //aux doubles

	double theta;		
	double sigma;
	double synp_eff;
	double tau;
	double input_det;
	double input_sto;
	double threshold;
	double Activity;

	/*Equations parameters*/
	double t=t_0_sim; 				//Initial time
	double x_0=0.0;				//Initial condition for m_ij
	int temp_steps= floor((t_f_sim-t_0_sim)/dt);  	// # of steps (time).
	

	/*Initial condition for m_{i->j} is zero*/
	int ***ocu=NULL;
	ocu=(int***)malloc(N_BGconnection*sizeof(int**));
	for(i=0;i<N_BGconnection;i++){
    	N_aux=Network[i].last->neuron_number;
    	ocu[i]=(int**)malloc(N_aux*sizeof(int*));
    	for(j=0;j<N_aux;j++){
        	Kab=Numbers_Connection[i][j];
        	ocu[i][j]=(int*)malloc(Kab*sizeof(int));
        	for(k=0;k<Kab;k++){
            	(Network[i].cond)[j][k].sinaptic_curr=x_0;
            	ocu[i][j][k]=dim_delay[i];}
    	}
	}


	/*Files - Results*/
	//FILE *SPIKES; 					//Save (neuron,time) of spike
	FILE *CTX, *STN,*GPI;//,*THA,*CTX,*STR; 	//Save current in BG.
	FILE *HCTX,*HSTR;  				//Verification External Input
	FILE *HDBS;						//Verification DBS
	//FILE *Output_A_B;

	//SPIKES=fopen("SPIKES.dat","w");

	STN=fopen("STN.dat","w");
	GPI=fopen("GPI.dat","w");
	//THA=fopen("THA.dat","w");
	CTX=fopen("CTX.dat","w");
	//STR=fopen("STR.dat","w");

	HCTX=fopen("H_CTX.dat","w");
	HSTR=fopen("H_STR.dat","w");
	HDBS=fopen("H_DBS.dat","w");

	//Output_A_B=fopen("Output_A_B.dat","w");

	/*DBS parameters*/
	int BG_dbs_index= floor(Electr_parameter[0]); //index_stimulation

	/*Random pulses*/
	//double frec_m=1.0/DBS_parameter[1];		//mean frequency
	//double per=Variability; 				//sigma/f_med
	//double interval=(t_f_DBS-t_0_DBS); 		//ms
	//int n_pulse=floor((frec_m*interval)*(1-pow(per,2)));
	//double *inicios=(double*)malloc(n_pulse*sizeof(double));

	//srand(time(NULL));
	//inicios[0]=t_0_DBS;
	
	//for(i=1;i<n_pulse;i++){
 	//   inicios[i]=inicios[i-1]+1.0/random_frequency(frec_m,per);}


	/*Temporal solution of eqs dif*/

	for(i=0;i<temp_steps;i++){
        fprintf(STN,"%f \t",t);
        fprintf(GPI,"%f \t",t);
        //fprintf(THA,"%f \t",t);
        fprintf(CTX,"%f \t",t);
        //fprintf(STR,"%f \t",t);
        fprintf(HCTX,"%f \t",t);
        fprintf(HSTR,"%f \t",t);
        fprintf(HDBS,"%f \t",t);
  		//fprintf(Output_A_B,"%f \t",t);

        /*----------------------------------------------------------------------------------------------------------*/
        /*Load external (deterministic and noise) input */
        for(j=0;j<Number_BG;j++){
                for(k=0;k<Population[j].neuron_number;k++){
                	theta=Population[j].circle[k].angle;
                    Population[j].circle[k].input_det = Population[j].circle[k].ext(ext_parameters[j],cos(theta),t);

                    sigma=Population[j].circle[k].parameters[2];
                    Population[j].circle[k].input_sto = external_noise(0,sigma);}
        }


        /*Total inputs*/
        for(j=0;j<N_BGconnection;j++){
                N_aux=Network[j].last->neuron_number;
                synp_eff=Network[j].parameters[0];

                for(k=0;k<N_aux;k++){
                        for(s=0;s<Numbers_Connection[j][k];s++){
                            borrarinicio(data_delayed[j][k][s],&(ocu[j][k][s]),&value);
                            ((Network[j].cond)[k][s].next)->input_det = ((Network[j].cond)[k][s].next)->input_det + (synp_eff*value);}
                }
        }


        /*DBS Stimulation*/
        for(k=0;k<Population[BG_dbs_index-1].neuron_number;k++){
        	theta=Population[BG_dbs_index-1].circle[k].angle;
            Population[BG_dbs_index-1].circle[k].input_det=Population[BG_dbs_index-1].circle[k].input_det+spatial_electrode(Electr_parameter[1],Electr_parameter[2],Electr_parameter[3],theta)*DBS_Periodic(t,DBS_parameter[0],DBS_parameter[1],DBS_parameter[2],DBS_parameter[3]);}
            //Population[BG_dbs_index-1].circle[k].input_det=Population[BG_dbs_index-1].circle[k].input_det+spatial_electrode(Electr_parameter[1],Electr_parameter[2],Electr_parameter[3],theta)*DBS_Random(t,DBS_parameter[0],DBS_parameter[2],DBS_parameter[3],inicios,n_pulse);}

		/*----------------------------------------------------------------------------------------------------------*/
        
        /*Verification external iputs*/
        fprintf(HCTX,"%f \n",Population[1].circle[1400].ext(ext_parameters[1],cos(Population[1].circle[1400].angle),t));
        fprintf(HSTR,"%f \n",Population[3].circle[1400].ext(ext_parameters[3],cos(Population[3].circle[1400].angle),t));
		fprintf(HDBS,"%f \n",DBS_Periodic(t,DBS_parameter[0],DBS_parameter[1],DBS_parameter[2],DBS_parameter[3]));
        //fprintf(HDBS,"%f \n",DBS_Random(t,DBS_parameter[0],DBS_parameter[2],DBS_parameter[3],inicios,n_pulse));


        /*Synaptic current*/
        /*for(k=0;k<Population[0].neuron_number;k++){
                fprintf(THA,"%lf \t",Population[0].circle[k].input_det+Population[0].circle[k].input_sto);}

		fprintf(THA,"\n");*/

        for(k=0;k<Population[1].neuron_number;k++){
                fprintf(CTX,"%lf \t",Population[1].circle[k].input_det+Population[1].circle[k].input_sto);}

		fprintf(CTX,"\n");
        
        for(k=0;k<Population[2].neuron_number;k++){
                fprintf(STN,"%lf \t",Population[2].circle[k].input_det+Population[2].circle[k].input_sto);}
		
		fprintf(STN,"\n");

        /*for(k=0;k<Population[3].neuron_number;k++){
                fprintf(STR,"%lf \t",Population[3].circle[k].input_det+Population[3].circle[k].input_sto);}
		
		fprintf(STR,"\n");*/

        for(k=0;k<Population[4].neuron_number;k++){
                fprintf(GPI,"%lf \t",Population[4].circle[k].input_det+Population[4].circle[k].input_sto);}
        
        fprintf(GPI,"\n");
        
        

        /*----------------------------------------------------------------------------------------------------------*/

        /*Generation of spikes*/
        /*for(k=0;k<Population[BG_Spike-1].neuron_number;k++){
               value=rand()/(double)RAND_MAX;
               input_det=Population[BG_Spike-1].circle[k].input_det;
               input_sto=Population[BG_Spike-1].circle[k].input_sto;
               threshold=Population[BG_Spike-1].circle[k].parameters[0];
               
               if(value< dt*(Activation_function(input_det,input_sto,threshold,&aux)+aux)){
                   fprintf(SPIKES,"%d \t %lf \n",k,t);}
        }

        /*----------------------------------------------------------------------------------------------------------*/

        /*VARIABLE M_Alpha->Beta (Output)*/
        //for (j=0;j<C;j++){
        //    fprintf(Output_A_B,"%f \t",(Red[j].cond)[0][0].M);}

        //fprintf(Output_A_B,"\n");

        /*----------------------------------------------------------------------------------------------------------*/

        /*Solving system*/
        for (j=0;j<N_BGconnection;j++){
                N_aux=Network[j].last->neuron_number;
                for(k=0;k<N_aux;k++){
                        for(s=0;s<Numbers_Connection[j][k];s++){
                        	input_det=(Network[j].cond)[k][s].last->input_det;
                        	input_sto=(Network[j].cond)[k][s].last->input_sto;
                            threshold=(Network[j].cond)[k][s].last->parameters[0];
                            tau=(Network[j].parameters[1]);
                            Activity=Activation_function(input_det,input_sto,threshold,&aux);

                            (Network[j].cond)[k][s].sinaptic_curr = (Network[j].cond)[k][s].sinaptic_curr + dt * (Activity -(Network[j].cond)[k][s].sinaptic_curr)/ tau + (sqrt(dt)*aux)/tau;
                            insertarfinal(data_delayed[j][k][s],&(ocu[j][k][s]),(Network[j].cond)[k][s].sinaptic_curr,dim_delay[j]);}
                }
        }

		t=t+dt;
	}

	//fclose(SPIKES);
	fclose(STN);
	fclose(GPI);
	//fclose(THA);
	fclose(CTX);
	//fclose(STR);
	fclose(HCTX);
	fclose(HSTR);
	fclose(HDBS);
	//fclose(Output_A_B);
}




int main(){
//int main(int argc, char *argv[]){ 
//argc = ["main", H_ctx, H_str,amplitude_DBS,period_DBS,width_DBS,slope_DBS,variability (option)]
	
	srand(time(NULL));
	char *endptr;
	
	int i,j,k;						//indexs
	int aux;

	/*-------------------BUILD NETWORK-------------------------------*/
	/*---------------------------------------------------------------*/	
	BGan *Population=NULL;     		//Populations of basal ganglia
	conductor *Network=NULL;		//Connections between neurons
	int **Numbers_Connection=NULL;	//Number of connection for each neuron (in each BG)
	int Number_BG=5;

	printf("1-Iniciando. \n");

	int N_BGconnection=Network_construction(&Population,Number_BG,&Network,&Numbers_Connection); //Building network - N_BGconnection: number of connections between BG's.

	printf("2-Red Construida. \n");

	/*--------------------DELAYED DATA-------------------------------*/
	/*---------------------------------------------------------------*/
	double ****data_delayed=(double ****)malloc(N_BGconnection*sizeof(double***)); //values of m_ij "t in(-delay,0)" [BG_index][neuron_index][connection_index][time] 	
	int dim_delay[N_BGconnection]; 												  //Dimension of delay vector for each BG connections
	

	for(i=0;i<N_BGconnection;i++){									  //for each connection i between BG's
    	dim_delay[i]=floor(Network[i].parameters[2]/dt); 			  //dim _ delay = delay / dt 
    	aux=Network[i].last->neuron_number;						      //aux:number of neurons in BG input 
    	data_delayed[i]=(double***)malloc(aux*sizeof(double**));	
    	
    	for(j=0;j<aux;j++){										//for each neuron j in this BG.
            data_delayed[i][j]=(double**)malloc(Numbers_Connection[i][j]*sizeof(double*)); 
            
            for(k=0;k<Numbers_Connection[i][j];k++){			//for each connection k of neuron j
                data_delayed[i][j][k]=(double*)calloc(dim_delay[i],sizeof(double));} //data delayed = vector of dim dim_delay[i]
            
    	}
	}

	printf("3-Vectores delays listos.\n");

	//- PARA VER DISTRIBUCION DE CONEXIONES -//
	//FILE *fp;
	//fp=fopen("Numbers_Connection-Gpi-Th.dat","w");
	//
	//for(i=0;i<2800;i++){
	//    fprintf(fp,"%d \n",Numbers_Connection[5][i]);
	//}
	//
	//fclose(fp);


	/*-----INPUTS EXTERNOS CTX-STR------*/
	double **ext_parameters=(double**)malloc(Number_BG*sizeof(double*)); 						//Numeration 1-Th; 2-Ctx; 3-STN; 4-Str; 5-GPi
	for (i=0;i<Number_BG;i++){																	//Mean spat, mod spat, mean temp, mod temp, period
		ext_parameters[i]=(double*)calloc(5,sizeof(double));
		ext_parameters[i][2]=1.0;
		ext_parameters[i][4]=1.0;
	}

	scanf("%lf \t %lf \t",&ext_parameters[1][0],&ext_parameters[3][1]);

	//ext_parameters[1][0]=strtod(argv[1],&endptr); //CTX
	//ext_parameters[3][1]=strtod(argv[2],&endptr); //Str
	
	for(j=0;j<Population[1].neuron_number;j++){
    	Change_InputExt(&Population,external_Ctx,2,j);}

	for(j=0;j<Population[3].neuron_number;j++){
    	Change_InputExt(&Population,external_Str,4,j);}

	printf("4-Se cambiaron las entradas.\n");


	/*----------------LFP PARAMETERS(geometric)-----------------------*/
	/*---------------------------------------------------------------*/

	FILE *Electr_LFP;
	Electr_LFP=fopen("ElectrodosLFP.dat","r");

	double **LFP_parameter=(double**)malloc(Number_BG*sizeof(double*));
	for(i=0;i<Number_BG;i++){										//for each BG: [index_BG,ang_center,threshold,volumen_integration] 
        LFP_parameter[i]=(double*)malloc(4*sizeof(double));
        fscanf(Electr_LFP,"%lf \t %lf \t %lf \t %lf \n",&LFP_parameter[i][0],&LFP_parameter[i][1],&LFP_parameter[i][2],&LFP_parameter[i][3]);}

	fclose(Electr_LFP);

	/*-------------STIMULATOR PARAMETERS(geometric)------------------*/
	/*---------------------------------------------------------------*/

	FILE *Electr_Geometry;
	Electr_Geometry=fopen("ElectrodosEstim.dat","r");

	double *Electr_parameter=(double*)malloc(4*sizeof(double));		//[index_BG,ang_center,threshold,VTA]
	fscanf(Electr_Geometry,"%lf \t %lf \t %lf \t %lf \n",&Electr_parameter[0],&Electr_parameter[1],&Electr_parameter[2],&Electr_parameter[3]); 

	fclose(Electr_Geometry);

	/*-------------STIMULATOR PARAMETERS(temporal)-------------------*/
	/*---------------------------------------------------------------*/

	double DBS_parameter[4]; 										//[amplitude,period,width,slope]
	scanf("%lf \t %lf \t %lf \t %lf",&DBS_parameter[0],&DBS_parameter[1],&DBS_parameter[2],&DBS_parameter[3]);

	//for(i=0;i<4;i++){
    //	DBS_parameter[i]=strtod(argv[i+3],&endptr);} 
	
	//double Variability=strtod(argv[7],&endptr);        			//variability in random DBS.

	printf("5-Electrodos configurados.\n");

	/*-----------------------SOLVING SYSTEM--------------------------*/
	/*---------------------------------------------------------------*/

	printf("6-Resolviendo sistema.\n");

	int BG_spike=3; 												//Numeration 1-Th; 2-Ctx; 3-STN; 4-Str; 5-GPi
	solve_equations(Population,Number_BG,Network,N_BGconnection,data_delayed,dim_delay,Numbers_Connection,BG_spike,ext_parameters,LFP_parameter,Electr_parameter,DBS_parameter);
	//solve_equations(Population,Number_BG,Network,N_BGconnection,data_delayed,dim_delay,Numbers_Connection,BG_spike,ext_parameters,LFP_parameter,Electr_parameter,DBS_parameter,Variability);

	printf("7-Fin\n");
	return 0;
}
