//////////////////////////////////////////////////
// file:        LDcodedMIMO.cpp                     //
// copywrite:   Philippe Bergeron-Burns         //
// description: Testing the 2 transmit antenna  //
//              LD and VBLAST                   //
//////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream.h>
#include <complex>
#include "matrix.cpp"
using namespace std;
using namespace math;
#define STD std
typedef complex<double> Complex;
typedef matrix<double> Matrix;
double pi = 3.1415926535;
#include "utilities.cpp"
#include "QR.cpp"
#include "sphere1.cpp"
#include "sphere2.cpp"

int main()
{
	///////////////  Define user constants //////////////////////////////////
	const int ACTUAL_NUM_RX = 2; //The number of transmit antennas          
	const int ACTUAL_NUM_TX = 2;  //The number of receive antennas          
	const int PAM = 2;     //The constellation.  supported: 2,4,8,16      
	const int TRIALS = 100; //The number of trials at each SNR.                   
	const int MIN_ERROR = 1; //the minimum of symbol errors at any SNR	

	//the two front end parameters:
	const double pure_Ph1 = 1.8;
	const double pure_Ph2 = 2.2;
	const double pure_Ph3 = 2.6;
	const double pure_Ph4 = 3.0;
	const double pure_Ph5 = 3.4;
	const double SPC_NOISE = 2.4; //The SPC-ML scheme variable U.

	for(int SNRdB = 12;SNRdB<=12;SNRdB +=2) {

	srand( (unsigned)time( NULL ) );  //seed the random number generator
	clock_t   start, finish;
	start = clock();                   //start the clock
	const int NUM_RX = 4 * ACTUAL_NUM_RX;
	const int NUM_TX = 4 * ACTUAL_NUM_TX;
	double SNR =  pow(10, (double)SNRdB / 10 );
	double num_bitsPERsymbol = 2*log((double)PAM)/log((double)2);
	double scale,snr; int i;int SCALE = 0;
	getScale(scale,SCALE,snr,SNR,ACTUAL_NUM_TX,PAM);
	double snr_db = 10 * log10( snr );
	int qam = PAM * PAM;
	printf("\n%d Trials. dB %f, scale %f, QAM-%d, %d --> %d\n",
		TRIALS,snr_db,scale,qam,ACTUAL_NUM_TX,ACTUAL_NUM_RX);	

	int xFixed[NUM_TX+1];  //SPC record keeping
	for(i=0;i<=NUM_TX;i++) xFixed[i] = 0;
	int numFixed = 0;      //number of total symbols the PC-ML method verified

	int prob_SPC = 0;
	int errors_ld = 0;	    
	int errors_bits_ld = 0; 
	int errors_v  = 0;	    
	int errors_bits_v  = 0; 
	int errors_v_v  = 0;	
	int errors_bits_v_v  = 0;
	int errors_ld_v  = 0;	 
	int errors_bits_ld_v  = 0;
	int errors_se  = 0;	      
	int errors_bits_se  = 0;  
	int errors_eq  = 0;	      
	int errors_bits_eq  = 0;  
	int errors_SPC_SD  = 0;	  
	int errors_bits_SPC_SD  = 0; 
	int errors_SPC_SE  = 0;	     
	int errors_bits_SPC_SE  = 0; 
	int loop_total_v = 0;
	int loop_total_ld = 0;
	int se_loop_total = 0;
	int loop_total_SPC_SD = 0;
	int loop_total_SPC_SE = 0;
	int percent1 = 0; int percent2 = 0; int percent3 = 0;
	int percent4 = 0; int percent5 = 0;
	int pure_errors1 = 0;
	int pure_errors2 = 0;
	int pure_errors3 = 0;
	int pure_errors4 = 0;
	int pure_errors5 = 0;
	int pure_errors_bits1 = 0;
	int pure_errors_bits2 = 0;
	int pure_errors_bits3 = 0;
	int pure_errors_bits4 = 0;
	int pure_errors_bits5 = 0;
	int pureSPC_SD_errors1 = 0;
	int pureSPC_SD_errors2 = 0;
	int pureSPC_SD_errors3 = 0;
	int pureSPC_SD_errors4 = 0;
	int pureSPC_SD_errors5 = 0;
	int pureSPC_SD_errors_bits1 = 0;
	int pureSPC_SD_errors_bits2 = 0;
	int pureSPC_SD_errors_bits3 = 0;
	int pureSPC_SD_errors_bits4 = 0;
	int pureSPC_SD_errors_bits5 = 0;
	unsigned long loop_total_Ph1_SD = 0;
	unsigned long loop_total_Ph2_SD = 0;
	unsigned long loop_total_Ph3_SD = 0;
	unsigned long loop_total_Ph4_SD = 0;
	unsigned long loop_total_Ph5_SD = 0;
	unsigned long loop_total_Ph1_SE = 0;
	unsigned long loop_total_Ph2_SE = 0;
	unsigned long loop_total_Ph3_SE = 0;
	unsigned long loop_total_Ph4_SE = 0;
	unsigned long loop_total_Ph5_SE = 0;
	unsigned long loop_total_Ph1_SPC_SD = 0;
	unsigned long loop_total_Ph2_SPC_SD = 0;
	unsigned long loop_total_Ph3_SPC_SD = 0;
	unsigned long loop_total_Ph4_SPC_SD = 0;
	unsigned long loop_total_Ph5_SPC_SD = 0;
	unsigned long loop_total_Ph1_SPC_SE = 0;
	unsigned long loop_total_Ph2_SPC_SE = 0;
	unsigned long loop_total_Ph3_SPC_SE = 0;
	unsigned long loop_total_Ph4_SPC_SE = 0;
	unsigned long loop_total_Ph5_SPC_SE = 0;
	Matrix eye_RX = Make_eye(NUM_RX, 1);	

	//The Monte-Carlo Loop
	int p = -1;
	while(errors_ld < MIN_ERROR || p < TRIALS)
	{
		p++;

		int numDim = 0;
		int num_ML = 0;			// number of ML symbol errors in the frame
		int num_eq = 0;			// number of eq symbol errors in the frame
		int num_SPC_SD = 0;
		int num_ML_bits = 0;	// number of ML bit errors in the frame
		int num_eq_bits = 0;	// number of eq bit errors in the frame
		int num_SPC_SD_bits = 0;


		Matrix H = makeChannel(ACTUAL_NUM_RX,ACTUAL_NUM_TX);
		Matrix H_ld = makeLD22_channel(H);
		Matrix H_v = makeV22_channel(H);

		Matrix ZF_ld = (!( (~H_ld) * H_ld )) * (~H_ld);
		Matrix ZF_v  = (!( (~H_v ) * H_v  )) * (~H_v );
		Matrix PH_ld = H_ld * ZF_ld;	

		Matrix bits = getInputs(NUM_TX,PAM);
		Matrix u(NUM_TX,1);
		for (i = 0; i < NUM_TX;i++)
			u(i,0) = bits(i,0) * scale;
		Matrix n = makeNoise(NUM_RX);

		Matrix y_ld = H_ld * u + n;
		Matrix y_v  = H_v  * u + n;

		Matrix Equ_ld = ZF_ld * y_ld;
		Matrix Equalized_ld(NUM_RX,1);
		for(i=0;i<NUM_TX;i++)
			Equalized_ld(i,0) = Equ_ld(i,0)/scale;
		Matrix Eq_bits_ld = get_Eq_bits(NUM_TX,Equalized_ld,PAM);
		Matrix difference_ld = Equalized_ld - Eq_bits_ld;
		int tester = 0;  //will count how many times the SD loops
		Matrix best_bits_ld = sphere_decoder(H_ld, NUM_TX, NUM_RX, 
			Equalized_ld, difference_ld, PAM, tester);
		loop_total_ld += tester;

		int se_tester = 0;
		Matrix se_best_bits = se_sphere_decoder(H_ld, NUM_TX, NUM_RX, 
			Equalized_ld, difference_ld, PAM, se_tester);
		se_loop_total += se_tester;

		Matrix Equ_v  = ZF_v  * y_v;
		Matrix Equalized_v(NUM_RX,1);
		for(i=0;i<NUM_TX;i++)
			Equalized_v(i,0) = Equ_v(i,0)/scale;
		Matrix Eq_bits_v = get_Eq_bits(NUM_TX,Equalized_v,PAM);
		Matrix difference_v = Equalized_v - Eq_bits_v;
		tester = 0;
		Matrix best_bits_v = sphere_decoder(H_v, NUM_TX, NUM_RX, 
			Equalized_v, difference_v, PAM, tester);
		loop_total_v += tester;

		Matrix VBLAST_OUT = vblast_decoder(H_v, NUM_TX, NUM_RX, y_v, scale, PAM);
		Matrix VBLAST_OUT_ld = vblast_decoder(H_ld, NUM_TX, NUM_RX, y_ld, scale, PAM);

		Matrix scaled_eq_bits(NUM_TX,1);
		for(i=0;i<NUM_TX;i++)
			scaled_eq_bits(i,0) = scale * Eq_bits_ld(i,0);
		Matrix pseudo_sent(NUM_RX,1);
		pseudo_sent = H_ld * scaled_eq_bits;
		Matrix w(NUM_RX,1);
		for(i=0;i<NUM_RX;i++)
		w(i,0) = y_ld(i,0) - pseudo_sent(i,0);
		Matrix SPC_u(NUM_TX,1);
		for(i = 0;i< NUM_TX; i++) SPC_u(i,0) = 0;

		for(int dim = 0;dim<NUM_TX;dim++) /*FOR EACH DIM*/
		{
			Matrix G(NUM_RX,1);
			for(i=0; i < NUM_RX;i++)
				G(i,0) = H_ld(i,dim);
			Matrix S(NUM_RX,NUM_TX-1);
			int reach = 0;
			for(int t = 0; t< NUM_TX; t++)
			{
				if(t == dim)
					reach = 1;
				else
				{
               		for(i=0;i<NUM_RX;i++)
              			S(i, t - reach) = H_ld(i,t);
				}
			}
			Matrix PS = S * (!( (~S) * S )) * (~S);		
			Matrix PS_n = eye_RX - PS;
			Matrix V = PS_n * G;
			Matrix PV = V * (!( (~V) * V )) * (~V);
			Matrix PVw = PV * w;
			double Vsize = eSize(V,NUM_RX);
			double PVwsize = eSize(PVw,NUM_RX);

			double theta = 2 * scale * Vsize - SPC_NOISE;
			if(PVwsize < theta)
			{
				SPC_u(dim,0) = Eq_bits_ld(dim,0);
				if(Eq_bits_ld(dim,0) != best_bits_ld(dim,0))   prob_SPC++;
				numFixed++;
				numDim++;
			}
		}  /*end of of FOR EACH DIM */
		Matrix SPC_SD = SPC_u;
		Matrix SPC_SE = SPC_u;
		Matrix y_SPC(NUM_RX, 1);
		y_SPC = y_ld;
		for(i = 0;i< NUM_TX; i++)
		{
			for(int j = 0; j< NUM_RX; j++)
				y_SPC(j,0) = y_SPC(j,0) - scale * H_ld(j,i) * SPC_u(i,0);
		}		
		int dim_left = NUM_TX - numDim;
		Matrix H_SPC(NUM_RX,dim_left);
		int tester_SPC_SD = 0; int tester_SPC_SE = 0;
		if(dim_left)
		{
			int upto = 0;
			for(i = 0;i< NUM_TX; i++)
			{
				if(SPC_u(i,0) == 0)
				{
					for(int j = 0;j < NUM_RX; j++)
						H_SPC(j,upto) = H_ld(j,i);
					upto++;
				}
			}
			Matrix ZF_SPC = (!( (~H_SPC) * H_SPC )) * (~H_SPC);
			Matrix Equ_SPC = ZF_SPC * y_SPC;
			Matrix Equalized_SPC(dim_left,1);
			Matrix Eq_bits_SPC(dim_left,1);
			for(i=0;i<dim_left;i++)
				Equalized_SPC(i,0) = Equ_SPC(i,0)/scale;
			Eq_bits_SPC = get_Eq_bits(dim_left,Equalized_SPC,PAM);
			Matrix difference_SPC = Equalized_SPC - Eq_bits_SPC;
			tester_SPC_SD = 0;
			Matrix best_bits_SPC_SD = sphere_decoder(H_SPC, dim_left, NUM_RX, 
				Equalized_SPC, difference_SPC, PAM, tester_SPC_SD);
			loop_total_SPC_SD += tester_SPC_SD;
			tester_SPC_SE = 0;
			Matrix best_bits_SPC_SE = se_sphere_decoder(H_SPC, dim_left, NUM_RX, 
				Equalized_SPC, difference_SPC, PAM, tester_SPC_SE);
			loop_total_SPC_SE += tester_SPC_SE;
	
			int SD_index = 0;
			int SE_index = 0;
			int noSD_index = 0;
			for(i = 0; i< NUM_TX; i++)
			{
				if(SPC_SD(i,0) == 0)
					SPC_SD(i,0) = best_bits_SPC_SD(SD_index++,0);
				if(SPC_SE(i,0) == 0)
					SPC_SE(i,0) = best_bits_SPC_SE(SE_index++,0);
			}
		}  //end of if(dim_left) for SPC decoder


		Matrix diff_bits_ld  = bits - best_bits_ld;
		Matrix diff_bits_v   = bits - best_bits_v ;
		Matrix diff_bits_v_v = bits - VBLAST_OUT;
		Matrix diff_bits_ld_v = bits - VBLAST_OUT_ld;
		Matrix diff_bits_se  = bits - se_best_bits;
		Matrix diff_bits_eq  = bits - Eq_bits_ld;
		Matrix diff_SPC_SD = bits - SPC_SD;
		Matrix diff_SPC_SE = bits - SPC_SE;
		int temp_bits = 0;
		for(i=0;i<NUM_TX;i++)  
		{
			if(diff_bits_ld(i,0) != 0)
			{
				temp_bits = count_bits(PAM,diff_bits_ld(i,0),best_bits_ld(i,0),bits(i,0));
				errors_bits_ld += temp_bits;
				errors_ld++;
				num_ML++;
				num_ML_bits += temp_bits;
			}
			if(diff_bits_v(i,0) != 0)
			{
				temp_bits = count_bits(PAM,diff_bits_v(i,0),best_bits_v(i,0),bits(i,0));
				errors_bits_v += temp_bits;
				errors_v++;
			}
			if(diff_bits_v_v(i,0) != 0)
			{
				temp_bits = count_bits(PAM,diff_bits_v_v(i,0),VBLAST_OUT(i,0),bits(i,0));
				errors_bits_v_v += temp_bits;
				errors_v_v++;
			}
			if(diff_bits_ld_v(i,0) != 0)
			{
				temp_bits = count_bits(PAM,diff_bits_ld_v(i,0),VBLAST_OUT_ld(i,0),bits(i,0));
				errors_bits_ld_v += temp_bits;
				errors_ld_v++;
			}
			if(diff_bits_se(i,0) != 0)
			{
				temp_bits = count_bits(PAM,diff_bits_se(i,0),se_best_bits(i,0),bits(i,0));
				errors_bits_se += temp_bits;
				errors_se++;
			}
			if(diff_bits_eq(i,0) != 0)
			{
				temp_bits = count_bits(PAM,diff_bits_eq(i,0),Eq_bits_ld(i,0),bits(i,0));
				errors_bits_eq += temp_bits;
				errors_eq++;
				num_eq++;
				num_eq_bits += temp_bits;				
			}
			if(diff_SPC_SD(i,0) != 0)
			{
				temp_bits = count_bits(PAM,diff_SPC_SD(i,0),SPC_SD(i,0),bits(i,0));
				errors_bits_SPC_SD += temp_bits;
				errors_SPC_SD++;
				num_SPC_SD++;
				num_SPC_SD_bits += temp_bits;
			}
			if(diff_SPC_SE(i,0) != 0)
			{
				temp_bits = count_bits(PAM,diff_SPC_SE(i,0),SPC_SE(i,0),bits(i,0));
				errors_bits_SPC_SE += temp_bits;
				errors_SPC_SE++;
			}
		}	
		for(i=0;i<4;i++)
		{
			if(diff_bits_ld(i,0) != 0 && diff_bits_ld(i+4,0) != 0)
			{	errors_ld--; num_ML--; }
			if(diff_bits_v(i,0) != 0 && diff_bits_v(i+4,0) != 0)
				errors_v--;
			if(diff_bits_v_v(i,0) != 0 && diff_bits_v_v(i+4,0) != 0)
				errors_v_v--;
			if(diff_bits_ld_v(i,0) != 0 && diff_bits_ld_v(i+4,0) != 0)
				errors_ld_v--;
			if(diff_bits_se(i,0) != 0 && diff_bits_se(i+4,0) != 0)
				errors_se--;
			if(diff_bits_eq(i,0) != 0 && diff_bits_eq(i+4,0) != 0)
			{	errors_eq--; num_eq--; }
			if(diff_SPC_SD(i,0) != 0 && diff_SPC_SD(i+4,0) != 0)
			{	errors_SPC_SD--; num_SPC_SD--; }
			if(diff_SPC_SE(i,0) != 0 && diff_SPC_SE(i+4,0) != 0)
				errors_SPC_SE--;

		}

		xFixed[numDim] = xFixed[numDim] + 1;

		Matrix w_pure = PH_ld * w;
		double size_pure = eSize(w_pure,NUM_RX);

		if(size_pure > pure_Ph1)
		{
			percent1++;
			pure_errors1 += num_ML;
			pure_errors_bits1 += num_ML_bits;
			pureSPC_SD_errors1 += num_SPC_SD;
			pureSPC_SD_errors_bits1 += num_SPC_SD_bits;
			loop_total_Ph1_SD += tester;
			loop_total_Ph1_SE += se_tester;
			loop_total_Ph1_SPC_SD += tester_SPC_SD;
			loop_total_Ph1_SPC_SE += tester_SPC_SE;
		}
		else
		{
			pure_errors1 += num_eq;
			pure_errors_bits1 += num_eq_bits;
			pureSPC_SD_errors1 += num_eq;
			pureSPC_SD_errors_bits1 += num_eq_bits;
		}
		if(size_pure > pure_Ph2)
		{
			percent2++;
			pure_errors2 += num_ML;
			pure_errors_bits2 += num_ML_bits;
			loop_total_Ph2_SD += tester;
			loop_total_Ph2_SE += se_tester;
			pureSPC_SD_errors2 += num_SPC_SD;
			pureSPC_SD_errors_bits2 += num_SPC_SD_bits;
			loop_total_Ph2_SPC_SD += tester_SPC_SD;
			loop_total_Ph2_SPC_SE += tester_SPC_SE;
		}
		else
		{
			pure_errors2 += num_eq;
			pure_errors_bits2 += num_eq_bits;
			pureSPC_SD_errors2 += num_eq;
			pureSPC_SD_errors_bits2 += num_eq_bits;
		}
		if(size_pure > pure_Ph3)
		{
			percent3++;
			pure_errors3 += num_ML;
			pure_errors_bits3 += num_ML_bits;
			loop_total_Ph3_SD += tester;
			loop_total_Ph3_SE += se_tester;
			pureSPC_SD_errors3 += num_SPC_SD;
			pureSPC_SD_errors_bits3 += num_SPC_SD_bits;
			loop_total_Ph3_SPC_SD += tester_SPC_SD;
			loop_total_Ph3_SPC_SE += tester_SPC_SE;
		}
		else
		{
			pure_errors3 += num_eq;
			pure_errors_bits3 += num_eq_bits;
			pureSPC_SD_errors3 += num_eq;
			pureSPC_SD_errors_bits3 += num_eq_bits;
		}
		if(size_pure > pure_Ph4)
		{
			percent4++;
			pure_errors4 += num_ML;
			pure_errors_bits4 += num_ML_bits;
			loop_total_Ph4_SD += tester;
			loop_total_Ph4_SE += se_tester;
			pureSPC_SD_errors4 += num_SPC_SD;
			pureSPC_SD_errors_bits4 += num_SPC_SD_bits;
			loop_total_Ph4_SPC_SD += tester_SPC_SD;
			loop_total_Ph4_SPC_SE += tester_SPC_SE;
		}
		else
		{
			pure_errors4 += num_eq;
			pure_errors_bits4 += num_eq_bits;
			pureSPC_SD_errors4 += num_eq;
			pureSPC_SD_errors_bits4 += num_eq_bits;
		}
		if(size_pure > pure_Ph5)
		{
			percent5++;
			pure_errors5 += num_ML;
			pure_errors_bits5 += num_ML_bits;
			loop_total_Ph5_SD += tester;
			loop_total_Ph5_SE += se_tester;
			pureSPC_SD_errors5 += num_SPC_SD;
			pureSPC_SD_errors_bits5 += num_SPC_SD_bits;
			loop_total_Ph5_SPC_SD += tester_SPC_SD;
			loop_total_Ph5_SPC_SE += tester_SPC_SE;
		}
		else
		{
			pure_errors5 += num_eq;
			pure_errors_bits5 += num_eq_bits;
			pureSPC_SD_errors5 += num_eq;
			pureSPC_SD_errors_bits5 += num_eq_bits;
		}


	} //end of TRIALS for loop

	printf("   LD & SD       SER %d %f BER %d %f loops %f\n",errors_ld,
		(double)errors_ld/(p*2*ACTUAL_NUM_TX),errors_bits_ld,(double)errors_bits_ld/
		(p*4*num_bitsPERsymbol),(double)loop_total_ld/p);
	printf("   LD & SE       SER %d %f BER %d %f loops %f\n",errors_se,
		(double)errors_se/(p*2*ACTUAL_NUM_TX),errors_bits_se,(double)errors_bits_se/
		(p*4*num_bitsPERsymbol),(double)se_loop_total/p);
	printf("   LD & SPC-SD   SER %d %f BER %d %f loops %f\n",errors_SPC_SD,
		(double)errors_SPC_SD/(p*2*ACTUAL_NUM_TX),errors_bits_SPC_SD,
		(double)errors_bits_SPC_SD/(p*4*num_bitsPERsymbol),(double)loop_total_SPC_SD/p);
	printf("   LD & SPC-SE   SER %d %f BER %d %f loops %f\n",errors_SPC_SE,
		(double)errors_SPC_SE/(p*2*ACTUAL_NUM_TX),errors_bits_SPC_SE,
		(double)errors_bits_SPC_SE/(p*4*num_bitsPERsymbol),(double)loop_total_SPC_SE/p);
	printf("   LD & VBLAST   SER %d %f BER %d %f\n",errors_ld_v,
		(double)errors_ld_v/(p*2*ACTUAL_NUM_TX),errors_bits_ld_v,
		(double)errors_bits_ld_v/(p*4*num_bitsPERsymbol));
	printf("   LD & ZF       SER %d %f BER %d %f\n",errors_eq,
		(double)errors_eq/(p*2*ACTUAL_NUM_TX),errors_bits_eq,
		(double)errors_bits_eq/(p*4*num_bitsPERsymbol));
	printf("   uncoded & SD  SER %d %f BER %d %f loops %f\n",errors_v,
		(double)errors_v/(p*2*ACTUAL_NUM_TX),errors_bits_v,
		(double)errors_bits_v/(p*4*num_bitsPERsymbol),(double)loop_total_v/TRIALS);
	printf("uncoded & VBLAST SER %d %f BER %d %f\n",errors_v_v,
		(double)errors_v_v/(p*2*ACTUAL_NUM_TX),errors_bits_v_v,
		(double)errors_bits_v_v/(p*4*num_bitsPERsymbol));

	printf("\nPh T=%f SER %f BER %f per %f loops %f SE loops %f\n", pure_Ph1, 
		(double)pure_errors1/(p*2*ACTUAL_NUM_TX),(double)pure_errors_bits1/
		(p*4*num_bitsPERsymbol),(double)(percent1)/p, (double)loop_total_Ph1_SD/p, 
		(double)loop_total_Ph1_SE/p);
	printf("Ph T=%f SER %f BER %f per %f loops %f SE loops %f\n", pure_Ph2, 
		(double)pure_errors2/(p*2*ACTUAL_NUM_TX),(double)pure_errors_bits2/
		(p*4*num_bitsPERsymbol),(double)(percent2)/p, (double)loop_total_Ph2_SD/p, 
		(double)loop_total_Ph2_SE/p);
	printf("Ph T=%f SER %f BER %f per %f loops %f SE loops %f\n", pure_Ph3, 
		(double)pure_errors3/(p*2*ACTUAL_NUM_TX),(double)pure_errors_bits3/
		(p*4*num_bitsPERsymbol),(double)(percent3)/p, (double)loop_total_Ph3_SD/p, 
		(double)loop_total_Ph3_SE/p);
	printf("Ph T=%f SER %f BER %f per %f loops %f SE loops %f\n", pure_Ph4, 
		(double)pure_errors4/(p*2*ACTUAL_NUM_TX),(double)pure_errors_bits4/
		(p*4*num_bitsPERsymbol),(double)(percent4)/p, (double)loop_total_Ph4_SD/p, 
		(double)loop_total_Ph4_SE/p);
	printf("Ph T=%f SER %f BER %f per %f loops %f SE loops %f\n", pure_Ph5, 
		(double)pure_errors5/(p*2*ACTUAL_NUM_TX),(double)pure_errors_bits5/
		(p*4*num_bitsPERsymbol),(double)(percent5)/p, (double)loop_total_Ph5_SD/p, 
		(double)loop_total_Ph5_SE/p);
	printf("Ph-SPC T=%f SER %f BER %f SD loops %f SE loops %f\n", pure_Ph1, 
		(double)pureSPC_SD_errors1/(p*2*ACTUAL_NUM_TX),(double)pureSPC_SD_errors_bits1/
		(p*4*num_bitsPERsymbol), (double)loop_total_Ph1_SPC_SD/p, 
		(double)loop_total_Ph1_SPC_SE/p);
	printf("Ph-SPC T=%f SER %f BER %f SD loops %f SE loops %f\n", pure_Ph2, 
		(double)pureSPC_SD_errors2/(p*2*ACTUAL_NUM_TX),(double)pureSPC_SD_errors_bits2/
		(p*4*num_bitsPERsymbol), (double)loop_total_Ph2_SPC_SD/p, 
		(double)loop_total_Ph2_SPC_SE/p);
	printf("Ph-SPC T=%f SER %f BER %f SD loops %f SE loops %f\n", pure_Ph3, 
		(double)pureSPC_SD_errors3/(p*2*ACTUAL_NUM_TX),(double)pureSPC_SD_errors_bits3/
		(p*4*num_bitsPERsymbol), (double)loop_total_Ph3_SPC_SD/p, 
		(double)loop_total_Ph3_SPC_SE/p);
	printf("Ph-SPC T=%f SER %f BER %f SD loops %f SE loops %f\n", pure_Ph4, 
		(double)pureSPC_SD_errors4/(p*2*ACTUAL_NUM_TX),(double)pureSPC_SD_errors_bits4/
		(p*4*num_bitsPERsymbol), (double)loop_total_Ph4_SPC_SD/p, 
		(double)loop_total_Ph4_SPC_SE/p);
	printf("Ph-SPC T=%f SER %f BER %f SD loops %f SE loops %f\n", pure_Ph5, 
		(double)pureSPC_SD_errors5/(p*2*ACTUAL_NUM_TX),(double)pureSPC_SD_errors_bits5/
		(p*4*num_bitsPERsymbol), (double)loop_total_Ph5_SPC_SD/p, 
		(double)loop_total_Ph5_SPC_SE/p);


	printf("SPC %f numFixed %d per %f prob %d",SPC_NOISE, numFixed,
		(double)numFixed/(p*NUM_TX),prob_SPC);
	for(i=0;i<=NUM_TX;i++) printf(" %d:%d",i,xFixed[i]);
	finish = clock();
	double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "\nTest took: %6.6f seconds\n", duration );
	printf("number of trials %d\n",p);
	} // for each SNR

	return 1;
}  //end of main()

