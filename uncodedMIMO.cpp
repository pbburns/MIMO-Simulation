/////////////////////////////////////////////////////
// file:        uncodedMIMO.cpp                    //
// copywrite:   Philippe Bergeron-Burns            //
// description: The uncoded simulation entry point //
//              and primary flow                   //
/////////////////////////////////////////////////////

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

//This is the program entry point
int main()
{
	///////////////  Define user constants /////////////////////
	const int ACTUAL_NUM_RX = 7; //The number of transmit antennas          
	const int ACTUAL_NUM_TX = 5;  //The number of receive antennas          
	const int PAM = 4;     //The constellation.  supported: 2,4,8,16      
	const int TRIALS = 100000; //The number of trials at each SNR.                   
	const int MIN_ERROR = 100; //the minimum of symbol errors at any SNR
	const double pure_Ph1 = 1.8;
	const double pure_Ph2 = 2.2;
	const double pure_Ph3 = 2.6;
	const double pure_Ph4 = 3.0;
	const double pure_Ph5 = 3.4;
	const double SPC_NOISE = 3.6; //The SPC-ML scheme variable U.
	const double VAR_UNCERT = 0.002;  //variance of the H estimation error
	///////////////   END of user constants ////////////////////////

	//set the range of signal-to-noise ratios to simulate at
	for(int SNRdB = 16;SNRdB<= 23;SNRdB +=1) {

	srand( (unsigned)time( NULL ) );  //seed the random number generator
	clock_t   start, finish;
	start = clock();                   //start the clock
	//set the number of effective real antennas
	const int NUM_RX = 2 * ACTUAL_NUM_RX;
	const int NUM_TX = 2 * ACTUAL_NUM_TX;
	double SNR =  pow(10, (double)SNRdB / 10 );
	//calculate the number of bits per symbol 
	// (each squared QAM symbol has two PAM symbols)
	double num_bitsPERsymbol = 2*log((double) PAM)/log((double) 2);
	int debug = 0;

	//initialize all the symbol and bit error counters
	int errors = 0;	    
	int errors_bits = 0;
	int se_errors = 0;	
	int se_errors_bits = 0;	
	int errors1 = 0;	
	int errors1_bits = 0;	
	int errors_v = 0;   
	int errors_v_bits = 0; 
	int errors_mmse = 0;   
	int eq_error = 0;  
	int eq_error_bits = 0; 
	int eq_error1 = 0;  
	int eq_ML_error = 0;
	int eq_not_ML = 0;  
	int problem = 0;    
	int pure_problem1 = 0;
	int pure_problem2 = 0;
	int pure_problem3 = 0;
	int pure_problem4 = 0;
	int pure_problem5 = 0;
	int miss_hits1 = 0; int miss_hits2 = 0;
	int hits1 = 0; int hits2 = 0;
	int pure_hits1 = 0;
	int pure_hits2 = 0;
	int pure_hits3 = 0;
	int pure_hits4 = 0;
	int pure_hits5 = 0;
	int pure_miss1 = 0;
	int pure_miss2 = 0;
	int pure_miss3 = 0;
	int pure_miss4 = 0;
	int pure_miss5 = 0;
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
	int pureSPC_SE_errors1 = 0;
	int pureSPC_SE_errors2 = 0;
	int pureSPC_SE_errors3 = 0;
	int pureSPC_SE_errors4 = 0;
	int pureSPC_SE_errors5 = 0;
	int prob_pcml2 = 0;   
	int rec_PCML_prob = 0;
	int errors_PC = 0;   
	int errors_SPC_equalized = 0;   
	int errors_SPC_SD = 0; 
	int errors_SPC_SD_bits = 0;
	int errors_SPC_SE = 0;   
	int errors_SPC_SE_bits = 0; 
	int i;  //looping variable
	int dim;

	int xFixed2[NUM_TX+1];  //SPC record keeping
	for(i=0;i<=NUM_TX;i++) xFixed2[i] = 0;

	int numFixed2 = 0;     //number of total symbols the SPC front-end verified
	int prob = 0;          //number of times the SPC was suboptimal for a symbol
	double scale,snr;
	int SCALE = 0;
	//set the scaling information symbol scaling factor based NUM_TX and SNR
	getScale(scale,SCALE,snr,SNR,ACTUAL_NUM_TX,PAM);
	double snr_db = 10 * log10( snr );
	int qam = PAM * PAM;
	printf("\n%d Trials. dB %f, scale %f, QAM-%d, %d --> %d\n",
		TRIALS,snr_db,scale,qam,ACTUAL_NUM_TX,ACTUAL_NUM_RX);	
	
	//initialize SD and loop counters to keep track of complexity
	unsigned long loop_total = 0;     
	unsigned long se_loop_total = 0;  
	unsigned long loop_total2 = 0;    
	unsigned long loop_total3 = 0;    
	unsigned long loop_total5 = 0;    
	unsigned long loop_total_SPC = 0; 
	unsigned long loop_total_noSPC = 0;
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
	unsigned long se2_loop_total = 0;
	Matrix eH(NUM_RX,NUM_TX);
	Matrix H(NUM_RX,NUM_TX);
	Matrix M(NUM_RX,NUM_TX);
	Matrix P(NUM_TX,NUM_TX);
	Matrix Pinv(NUM_TX,NUM_TX);
	Matrix ZF(NUM_TX,NUM_RX);
	Matrix MMSE(NUM_TX,NUM_RX);
	Matrix PH(NUM_RX,NUM_RX);
	double Gsize = 0;      
	double Gsize_vec[NUM_TX];
	double G1size_vec[NUM_TX];
	double Vsize_vec[NUM_TX];
	Matrix eye_RX = Make_eye(NUM_RX, 1);	
	Matrix eye_TX = Make_eye(NUM_TX, 1/scale);	
	
	//The Monte-Carlo Loop.
	int p = 0;
	while(errors < MIN_ERROR || p < TRIALS)
	{
		p++;
		//generate the realization of the random channel matrix
		eH = makeChannel(ACTUAL_NUM_RX,ACTUAL_NUM_TX);
	//	Print(eH);printf("\n");
		//add estimation errors to the channel
		H = addErrors( eH, ACTUAL_NUM_RX, ACTUAL_NUM_TX, VAR_UNCERT);
	//	Print(H);
		//calculate the ZF matrix
		ZF = (!( (~H) * H )) * (~H);
		//calculate the orthogonal projection matrix onto the <H> subspace
		PH = H * (!( (~H) * H )) * (~H);	
		Matrix base(NUM_RX,1);
		for(dim = 0;dim< NUM_TX; dim++)
		{
			for(i=0;i<NUM_RX;i++)
				base(i,0) = H(i,dim);
			Gsize_vec[dim] = eSize( base, NUM_RX );
			for(i=0;i<NUM_RX;i++)
				base(i,0) = M(i,dim);
			G1size_vec[dim] = eSize( base, NUM_RX );
		}
		int one_ML = 0;        // at least one SD error in the frame.
		int num_ML = 0;			// number of ML symbol errors in the frame
		int num_eq = 0;			// number of eq symbol errors in the frame
		int num_ML_bits = 0;	// number of ML bit errors in the frame
		int num_eq_bits = 0;	// number of eq bit errors in the frame
		int num_SPC_SD = 0;
		int num_SPC_SD_bits = 0;
		int num_SPC_SE = 0;
		int one_eq_ML = 0;    //any TX bit equalization different from ML bit
		int one_eq_ML1 = 0;    //any TX bit LLL equalization different from ML bit
		double size = 0;      //euclidean norm squared of the projection of w onto Hi
		double size_vec[NUM_TX];
		double size1 = 0;      //euclidean norm squared of the projection of w onto Hi
		Matrix PCMLbits(NUM_TX,1);
		Matrix dont_incl(NUM_TX,1);
		int num_ml_eq = 0;

		Matrix bits = getInputs(NUM_TX,PAM);
		if(debug)
		{
			printf("transmitted bits\n");
			Print( bits );
		}
		Matrix u(NUM_TX,1);
		for (i = 0; i < NUM_TX;i++)
			u(i,0) = bits(i,0) * scale;
		Matrix n = makeNoise(NUM_RX);
		double nsize = eSize(n,NUM_RX);
		if(debug) printf("the size of the noise is: %f\n", nsize);
		Matrix y = eH * u + n;

		Matrix Equ = ZF * y;
		Matrix Equalized(NUM_TX,1);
		for(i=0;i<NUM_TX;i++)
			Equalized(i,0) = Equ(i,0)/scale;
		//make a decision on each ZF symbols
		Matrix Eq_bits = get_Eq_bits(NUM_TX,Equalized,PAM);
		Matrix difference = Equalized - Eq_bits;
		int tester = 0;  //will count the SD loops
		Matrix best_bits = sphere_decoder(H, NUM_TX, NUM_RX, 
			               Equalized, difference, PAM, tester);
		loop_total += tester;
		int se_tester = 0;
		//now use the SE/SD
		Matrix se_best_bits = se_sphere_decoder(H, NUM_TX, NUM_RX, Equalized, 
			                  difference, PAM, se_tester);
		se_loop_total += se_tester;

		Matrix VBLAST_OUT = vblast_decoder(H, NUM_TX, NUM_RX, y, scale, PAM);
		
		//find the difference between frames so that the errors can be counted
		Matrix diff_bits = bits - best_bits;
		Matrix se_diff_bits = bits - se_best_bits;
		Matrix diff_eq = bits - Eq_bits;
		Matrix diff_eq_ML = best_bits - Eq_bits;
		Matrix diff_v = bits - VBLAST_OUT;
		for(i=0;i<NUM_TX;i++)  
		{
			int temp_bits = 0;
			if(diff_bits(i,0) != 0)
			{
				//convert from symbol errors to bit errors
				temp_bits = count_bits(PAM,diff_bits(i,0),best_bits(i,0),bits(i,0));
				errors_bits += temp_bits;
				one_ML = 1;
				num_ML++;
				num_ML_bits += temp_bits;
				errors++;
			}
			if(se_diff_bits(i,0) != 0)
			{
				//convert from symbol errors to bit errors
				temp_bits = count_bits(PAM,se_diff_bits(i,0),se_best_bits(i,0),bits(i,0));
				se_errors_bits += temp_bits;
				se_errors++;
			}
			if(diff_v(i,0) != 0)
			{
				//convert from symbol errors to bit errors
				temp_bits = count_bits(PAM,diff_v(i,0),VBLAST_OUT(i,0),bits(i,0));
				errors_v_bits += temp_bits;
				errors_v++;
			}
			if(diff_eq(i,0) != 0)
			{
				//convert from symbol errors to bit errors
				temp_bits = count_bits(PAM,diff_eq(i,0), bits(i,0), Eq_bits(i,0));
				eq_error_bits += temp_bits;
				num_eq++;
				num_eq_bits += temp_bits;				
				eq_error++;
			}
			if(diff_eq_ML(i,0) != 0)
			{
				one_eq_ML = 1;
				eq_ML_error++;
				num_ml_eq++;
			}
		}
		//take into account that for SER don't count both Im and Re errors
		for(i=0;i<ACTUAL_NUM_TX;i++)
		{
			if(diff_bits(i,0) != 0 && diff_bits(i+ACTUAL_NUM_TX,0) != 0)
			{
				errors--;
				num_ML--;
			}
			if(se_diff_bits(i,0) != 0 && se_diff_bits(i+ACTUAL_NUM_TX,0) != 0)
				se_errors--;
			if(diff_v(i,0) != 0 && diff_v(i+ACTUAL_NUM_TX,0) != 0)
				errors_v--;
			if(diff_eq(i,0) != 0 && diff_eq(i+ACTUAL_NUM_TX,0) != 0)
			{
				eq_error--;
				num_eq--;
			}
			if(diff_eq_ML(i,0) != 0 && diff_eq_ML(i+ACTUAL_NUM_TX,0) != 0)
				num_ml_eq--;
		}

		///////////////  Do the Subspace Matched Filtering ///////////

		Matrix scaled_eq_bits(NUM_TX,1);
		for(i=0;i<NUM_TX;i++)
			scaled_eq_bits(i,0) = scale * Eq_bits(i,0);

		Matrix pseudo_sent(NUM_RX,1);
		pseudo_sent = H * scaled_eq_bits;

		Matrix w(NUM_RX,1);
		for(i=0;i<NUM_RX;i++)
			w(i,0) = y(i,0) - pseudo_sent(i,0);

		int reach = 0;
		int numDim2 = 0;
		for(i = 0;i<NUM_TX; i++) dont_incl(i,0) = 0;
		Matrix SPC_u(NUM_TX,1);
		for(i = 0;i< NUM_TX; i++) SPC_u(i,0) = 0;
		for(dim = 0;dim<NUM_TX;dim++)   /**FOR EACH DIMENSION **/
		{
			Matrix G(NUM_RX,1);
			for(i=0;i<NUM_RX;i++)
				G(i,0) = H(i,dim);
			Gsize = Gsize_vec[dim]; 
		
			Matrix PG(NUM_RX,NUM_RX);
			Matrix w3(NUM_RX,1);
			PG = G * (!( (~G) * G )) * (~G);
			w3 = PG * w;
			size = eSize(w3,NUM_RX);
			size_vec[dim] = size;
	
			Matrix S(NUM_RX,NUM_TX-1);
			reach = 0;
			for(int t = 0; t< NUM_TX; t++)
			{
				if(t == dim)
					reach = 1;
				else
				{
               		for(i=0;i<NUM_RX;i++)
					{
              			S(i, t - reach) = H(i,t);
					}
				}
			}
			//calculate the projection matrix onto the S_i subspace
			Matrix PS = S * (!( (~S) * S )) * (~S);		
			//calculate the projection matrix onto the orthogonal S_i subspace
			Matrix PS_n = eye_RX - PS;
			Matrix V = PS_n * G;
			Matrix PV = V * (!( (~V) * V )) * (~V);
			Matrix PVw = PV * w;
			double Vsize = eSize(V,NUM_RX);
			Vsize_vec[dim] = Vsize;
			double PVwsize = eSize(PVw,NUM_RX);

			double topSize = 0;
			for(i=0; i< NUM_RX;i++)
				topSize = topSize + G(i,0) * V(i,0);
			double cos_theta = topSize/(Vsize * Gsize);
			double ang = acos(cos_theta);

			double sum_Hj = 0;
			double largest_Hj = 0;
			double sum_Mj = 0;
			double largest_Mj = 0;
			for(i=0;i<NUM_TX;i++)
			{
				if(i != dim)  
				{
					sum_Hj = sum_Hj + Gsize_vec[i];
					if(Gsize_vec[i] > largest_Hj) largest_Hj = Gsize_vec[i];
				}
			}

			//**************** The SPC-ML scheme *************************/

			double theta2 = 2 * scale * Vsize - SPC_NOISE;
			if(debug) printf("theta2 is %f\n",theta2);
			if(PVwsize < theta2)
			{
				SPC_u(dim,0) = Eq_bits(dim,0);
				dont_incl(dim,0) = Eq_bits(dim,0);
				numFixed2++;
				numDim2++;
				if(Eq_bits(dim,0) != best_bits(dim,0))   prob_pcml2++;
			}

			//*****************END of SPC-ML scheme*********************/
			
		}  /*end of FOR EACH LATTICE DIMENSION  ************/

		/*  beginning of SPC rest of the symbols detection TESTING ************/
		
		Matrix SPC_equalized = SPC_u;
		Matrix SPC_SD = SPC_u;
		Matrix SPC_SE = SPC_u;
		Matrix noSPC_SD = SPC_u;
		int number_left = 0;
		Matrix y_SPC(NUM_RX, 1);
		y_SPC = y;
		for(i = 0;i< NUM_TX; i++)
		{
			for(int j = 0; j< NUM_RX; j++)
				y_SPC(j,0) = y_SPC(j,0) - scale * H(j,i) * SPC_u(i,0);
		}
		
		int dim_left = NUM_TX - numDim2;
		Matrix H_SPC(NUM_RX,dim_left);
		Matrix Equalized_SPC(dim_left,1);
		Matrix Equalized_noSPC(dim_left,1);
		Matrix Equ_SPC(dim_left,1);
		Matrix Eq_bits_SPC(dim_left,1);
		Matrix Eq_bits_noSPC(dim_left,1);
		Matrix ZF_SPC(dim_left,NUM_RX);
		int tester_SPC = 0;
		int se2_tester = 0;
		if(dim_left)
		{
			int upto = 0;
			for(i = 0;i< NUM_TX; i++)
			{
				if(SPC_u(i,0) == 0)
				{
					for(int j = 0;j < NUM_RX; j++)
					{
						H_SPC(j,upto) = H(j,i);
						Equalized_noSPC(upto,0) = Equalized(i,0);
						Eq_bits_noSPC(upto,0) = Eq_bits(i,0);
					}
					upto++;
				}
			}
			//calculate the ZF matrix excluding the pre-detected columns
			ZF_SPC = (!( (~H_SPC) * H_SPC )) * (~H_SPC);
			Equ_SPC = ZF_SPC * y_SPC;
			for(i=0;i<dim_left;i++)
				Equalized_SPC(i,0) = Equ_SPC(i,0)/scale;
			Eq_bits_SPC = get_Eq_bits(dim_left,Equalized_SPC,PAM);
			Matrix difference_SPC = Equalized_SPC - Eq_bits_SPC;
			Matrix difference_noSPC = Equalized_noSPC - Eq_bits_noSPC;
			int tester_noSPC = 0;
			Matrix best_bits_SPC = sphere_decoder(H_SPC, dim_left, NUM_RX, 
				        Equalized_SPC, difference_SPC, PAM, tester_SPC);
			Matrix best_bits_noSPC = best_bits;
			loop_total_SPC += tester_SPC;
			loop_total_noSPC += tester_noSPC;

			Matrix se2_best_bits = se_sphere_decoder(H_SPC, dim_left, 
				       NUM_RX, Equalized_SPC, difference_SPC, PAM, se2_tester);
			se2_loop_total += se2_tester;
			
			int eq_index = 0;
			int SD_index = 0;
			int SE_index = 0;
			int noSD_index = 0;
			for(i = 0; i< NUM_TX; i++)
			{
				if(SPC_equalized(i,0) == 0)
					SPC_equalized(i,0) = Eq_bits_SPC(eq_index++,0);
				if(SPC_SD(i,0) == 0)
					SPC_SD(i,0) = best_bits_SPC(SD_index++,0);
				if(SPC_SE(i,0) == 0)
					SPC_SE(i,0) = se2_best_bits(SE_index++,0);
				if(noSPC_SD(i,0) == 0)
					noSPC_SD(i,0) = best_bits_noSPC(noSD_index++,0);
			}
			
		}

		/*  end of SPC rest ************************************/

		Matrix diff_SPC_equalized = bits - SPC_equalized;
		Matrix diff_SPC_SD = bits - SPC_SD;
		Matrix diff_SPC_SE = bits - SPC_SE;
		for(i=0;i<NUM_TX;i++)
		{
			int temp_bits = 0;
			if(diff_SPC_equalized(i,0) != 0) errors_SPC_equalized++;
			if(diff_SPC_SD(i,0) != 0)
			{
				temp_bits = count_bits(PAM,diff_SPC_SD(i,0),SPC_SD(i,0),bits(i,0));
				errors_SPC_SD_bits += temp_bits;
				errors_SPC_SD++;
				num_SPC_SD++;
				num_SPC_SD_bits += temp_bits;
			}
			if(diff_SPC_SE(i,0) != 0) 
			{
				temp_bits = count_bits(PAM,diff_SPC_SE(i,0),SPC_SE(i,0),bits(i,0));
				errors_SPC_SE_bits += temp_bits;
				errors_SPC_SE++;
				num_SPC_SE++;
			}
		}
		//take into account that for SER don't count both Im and Re errors
		for(i=0;i<ACTUAL_NUM_TX;i++)
		{
			if(diff_SPC_equalized(i,0) != 0 && diff_SPC_equalized(i+ACTUAL_NUM_TX,0) != 0)
				errors_SPC_equalized--;
			if(diff_SPC_SD(i,0) != 0 && diff_SPC_SD(i+ACTUAL_NUM_TX,0) != 0)
				errors_SPC_SD--;
			if(diff_SPC_SE(i,0) != 0 && diff_SPC_SE(i+ACTUAL_NUM_TX,0) != 0)
				errors_SPC_SE--;
		}

		
		xFixed2[numDim2] = xFixed2[numDim2] + 1;

		//key step to SFC front-end: orthogonally project the 
		// decision feedback vector, w, onto the H subspace.
		Matrix w_pure = PH * w;
		double size_pure = eSize(w_pure,NUM_RX);

		//For this projection compare it to 5 different values of the T parameter
		//compare to 5 different values of T so as to obtain results for all 5 at
		//the same time.
		if(size_pure > pure_Ph1)
		{
			if(!one_eq_ML) pure_miss1++;
			else           pure_hits1++;				
			pure_errors1 += num_ML;
			pure_errors_bits1 += num_ML_bits;
			pureSPC_SD_errors1 += num_SPC_SD;
			pureSPC_SD_errors_bits1 += num_SPC_SD_bits;
			pureSPC_SE_errors1 += num_SPC_SE;
			loop_total_Ph1_SD += tester;
			loop_total_Ph1_SE += se_tester;
			loop_total_Ph1_SPC_SD += tester_SPC;
			loop_total_Ph1_SPC_SE += se2_tester;
		}
		else
		{
			if(one_eq_ML)
				pure_problem1 += num_ml_eq;
			pure_errors1 += num_eq;
			pure_errors_bits1 += num_eq_bits;
			pureSPC_SD_errors1 += num_eq;
			pureSPC_SD_errors_bits1 += num_eq_bits;
			pureSPC_SE_errors1 += num_eq;
		}
		if(size_pure > pure_Ph2)
		{
			if(!one_eq_ML) pure_miss2++;
			else           pure_hits2++;				
			pure_errors2 += num_ML;
			pure_errors_bits2 += num_ML_bits;
			loop_total_Ph2_SD += tester;
			loop_total_Ph2_SE += se_tester;
			pureSPC_SD_errors2 += num_SPC_SD;
			pureSPC_SD_errors_bits2 += num_SPC_SD_bits;
			pureSPC_SE_errors2 += num_SPC_SE;
			loop_total_Ph2_SPC_SD += tester_SPC;
			loop_total_Ph2_SPC_SE += se2_tester;
		}
		else
		{
			if(one_eq_ML)
				pure_problem2 += num_ml_eq;
			pure_errors2 += num_eq;
			pure_errors_bits2 += num_eq_bits;
			pureSPC_SD_errors2 += num_eq;
			pureSPC_SD_errors_bits2 += num_eq_bits;
			pureSPC_SE_errors2 += num_eq;
		}
		if(size_pure > pure_Ph3)
		{
			if(!one_eq_ML) pure_miss3++;
			else           pure_hits3++;				
			pure_errors3 += num_ML;
			pure_errors_bits3 += num_ML_bits;
			loop_total_Ph3_SD += tester;
			loop_total_Ph3_SE += se_tester;
			pureSPC_SD_errors3 += num_SPC_SD;
			pureSPC_SD_errors_bits3 += num_SPC_SD_bits;
			pureSPC_SE_errors3 += num_SPC_SE;
			loop_total_Ph3_SPC_SD += tester_SPC;
			loop_total_Ph3_SPC_SE += se2_tester;
		}
		else
		{
			if(one_eq_ML)
				pure_problem3 += num_ml_eq;
			pure_errors3 += num_eq;
			pure_errors_bits3 += num_eq_bits;
			pureSPC_SD_errors3 += num_eq;
			pureSPC_SD_errors_bits3 += num_eq_bits;
			pureSPC_SE_errors3 += num_eq;
		}
		if(size_pure > pure_Ph4)
		{
			if(!one_eq_ML) pure_miss4++;
			else           pure_hits4++;				
			pure_errors4 += num_ML;
			pure_errors_bits4 += num_ML_bits;
			loop_total_Ph4_SD += tester;
			loop_total_Ph4_SE += se_tester;
			pureSPC_SD_errors4 += num_SPC_SD;
			pureSPC_SD_errors_bits4 += num_SPC_SD_bits;
			pureSPC_SE_errors4 += num_SPC_SE;
			loop_total_Ph4_SPC_SD += tester_SPC;
			loop_total_Ph4_SPC_SE += se2_tester;
		}
		else
		{
			if(one_eq_ML)
				pure_problem4 += num_ml_eq;
			pure_errors4 += num_eq;
			pure_errors_bits4 += num_eq_bits;
			pureSPC_SD_errors4 += num_eq;
			pureSPC_SD_errors_bits4 += num_eq_bits;
			pureSPC_SE_errors4 += num_eq;
		}
		if(size_pure > pure_Ph5)
		{
			if(!one_eq_ML) pure_miss5++;
			else           pure_hits5++;
			pure_errors5 += num_ML;
			pure_errors_bits5 += num_ML_bits;
			loop_total_Ph5_SD += tester;
			loop_total_Ph5_SE += se_tester;
			pureSPC_SD_errors5 += num_SPC_SD;
			pureSPC_SD_errors_bits5 += num_SPC_SD_bits;
			pureSPC_SE_errors5 += num_SPC_SE;
			loop_total_Ph5_SPC_SD += tester_SPC;
			loop_total_Ph5_SPC_SE += se2_tester;
		}
		else
		{
			if(one_eq_ML)
				pure_problem5 += num_ml_eq;
			pure_errors5 += num_eq;
			pure_errors_bits5 += num_eq_bits;
			pureSPC_SD_errors5 += num_eq;
			pureSPC_SD_errors_bits5 += num_eq_bits;
			pureSPC_SE_errors5 += num_eq;
		}

		
	}  // for each TRIAL...

	//Print the error-performance and complexity results for the SNR being tested.
	//these results are shown for each of the detection methods.
	printf("          SD/sent SER %d %f BER %d %f loops %f\n",errors,
		(double)errors/(p*ACTUAL_NUM_TX),errors_bits,
		(double)errors_bits/(p*ACTUAL_NUM_TX*num_bitsPERsymbol),(double)loop_total/p);
	printf("          SE/sent SER %d %f BER %d %f loops %f\n",se_errors,
		(double)se_errors/(p*ACTUAL_NUM_TX),se_errors_bits,
		(double)se_errors_bits/(p*ACTUAL_NUM_TX*num_bitsPERsymbol),(double)se_loop_total/p);
	printf("      SPC-SD/sent SER %d %f BER %d %f loops %f\n",errors_SPC_SD,
		(double)errors_SPC_SD/(p*ACTUAL_NUM_TX),errors_SPC_SD_bits,
		(double)errors_SPC_SD_bits/(p*ACTUAL_NUM_TX*num_bitsPERsymbol),(double)loop_total_SPC/p);
	printf("      SPC-SE/sent SER %d %f BER %d %f loops %f\n",errors_SPC_SE,
		(double)errors_SPC_SE/(p*ACTUAL_NUM_TX),errors_SPC_SE_bits,
		(double)errors_SPC_SE_bits/(p*ACTUAL_NUM_TX*num_bitsPERsymbol),(double)se2_loop_total/p);
	printf("          ZF/sent SER %d %f BER %d %f\n",eq_error,
		(double)eq_error/(p*ACTUAL_NUM_TX),eq_error_bits,
		(double)eq_error_bits/(p*ACTUAL_NUM_TX*num_bitsPERsymbol));
	printf("   ZF VBLAST/sent SER %d %f BER %d %f\n",errors_v,
		(double)errors_v/(p*ACTUAL_NUM_TX),errors_v_bits,
		(double)errors_v_bits/(p*ACTUAL_NUM_TX*num_bitsPERsymbol));
	printf("SPC %f UNCERT %f numFixed %d per %f prob %d",SPC_NOISE, VAR_UNCERT,
		numFixed2,(double)numFixed2/(p*NUM_TX),prob_pcml2);
	for(i=0;i<=NUM_TX;i++) printf(" %d:%d",i,xFixed2[i]);
	printf("\nPh T=%f SER %f BER %f per %f loops %f SE loops %f\n", 
		pure_Ph1, (double)pure_errors1/(p*ACTUAL_NUM_TX),
		(double)pure_errors_bits1/(p*ACTUAL_NUM_TX*num_bitsPERsymbol),
		(double)(pure_hits1+pure_miss1)/p, (double)loop_total_Ph1_SD/p, 
		(double)loop_total_Ph1_SE/p);
	printf("Ph T=%f SER %f BER %f per %f loops %f SE loops %f\n",
		pure_Ph2, (double)pure_errors2/(p*ACTUAL_NUM_TX),
		(double)pure_errors_bits2/(p*ACTUAL_NUM_TX*num_bitsPERsymbol),
		(double)(pure_hits2+pure_miss2)/p, (double)loop_total_Ph2_SD/p, 
		(double)loop_total_Ph2_SE/p);
	printf("Ph T=%f SER %f BER %f per %f loops %f SE loops %f\n", pure_Ph3, 
		(double)pure_errors3/(p*ACTUAL_NUM_TX),(double)pure_errors_bits3/
		(p*ACTUAL_NUM_TX*num_bitsPERsymbol),(double)(pure_hits3+pure_miss3)/p, 
		(double)loop_total_Ph3_SD/p, (double)loop_total_Ph3_SE/p);
	printf("Ph T=%f SER %f BER %f per %f loops %f SE loops %f\n", pure_Ph4, 
		(double)pure_errors4/(p*ACTUAL_NUM_TX),(double)pure_errors_bits4/
		(p*ACTUAL_NUM_TX*num_bitsPERsymbol),(double)(pure_hits4+pure_miss4)/p, 
		(double)loop_total_Ph4_SD/p, (double)loop_total_Ph4_SE/p);
	printf("Ph T=%f SER %f BER %f per %f loops %f SE loops %f\n", pure_Ph5, 
		(double)pure_errors5/(p*ACTUAL_NUM_TX),(double)pure_errors_bits5/
		(p*ACTUAL_NUM_TX*num_bitsPERsymbol),(double)(pure_hits5+pure_miss5)/p, 
		(double)loop_total_Ph5_SD/p, (double)loop_total_Ph5_SE/p);
	printf("Ph-SPC T=%f SER %f BER %f SD loops %f SE loops %f\n", pure_Ph1, 
		(double)pureSPC_SD_errors1/(p*ACTUAL_NUM_TX),(double)pureSPC_SD_errors_bits1/
		(p*ACTUAL_NUM_TX*num_bitsPERsymbol), (double)loop_total_Ph1_SPC_SD/p, 
		(double)loop_total_Ph1_SPC_SE/p);
	printf("Ph-SPC T=%f SER %f BER %f SD loops %f SE loops %f\n", pure_Ph2, 
		(double)pureSPC_SD_errors2/(p*ACTUAL_NUM_TX),(double)pureSPC_SD_errors_bits2/
		(p*ACTUAL_NUM_TX*num_bitsPERsymbol), (double)loop_total_Ph2_SPC_SD/p, 
		(double)loop_total_Ph2_SPC_SE/p);
	printf("Ph-SPC T=%f SER %f BER %f SD loops %f SE loops %f\n", pure_Ph3, 
		(double)pureSPC_SD_errors3/(p*ACTUAL_NUM_TX),(double)pureSPC_SD_errors_bits3/
		(p*ACTUAL_NUM_TX*num_bitsPERsymbol), (double)loop_total_Ph3_SPC_SD/p, 
		(double)loop_total_Ph3_SPC_SE/p);
	printf("Ph-SPC T=%f SER %f BER %f SD loops %f SE loops %f\n", pure_Ph4, 
		(double)pureSPC_SD_errors4/(p*ACTUAL_NUM_TX),(double)pureSPC_SD_errors_bits4/
		(p*ACTUAL_NUM_TX*num_bitsPERsymbol), (double)loop_total_Ph4_SPC_SD/p, 
		(double)loop_total_Ph4_SPC_SE/p);
	printf("Ph-SPC T=%f SER %f BER %f SD loops %f SE loops %f\n", pure_Ph5, 
		(double)pureSPC_SD_errors5/(p*ACTUAL_NUM_TX),(double)pureSPC_SD_errors_bits5/
		(p*ACTUAL_NUM_TX*num_bitsPERsymbol), (double)loop_total_Ph5_SPC_SD/p, 
		(double)loop_total_Ph5_SPC_SE/p);
	finish = clock();
	double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "Test took: %6.6f seconds\n", duration );
	printf("number of trials %d\n",p);
	
	}  //SNRdB loop
	
	return 1;
}  // main()


