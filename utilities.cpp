
//Returns the norm of a column vector
double eSize(Matrix A, int NUM_RX)
{
	double s = 0;
	for(int i= 0; i< NUM_RX; i++)
		s = s + A(i,0) * A(i,0);
	s = sqrt(s);
	return s;
}

//Prints the elements of a matrix
void Print(Matrix M)
{
	for (int i=0; i < M.RowNo(); i++)
	{
      for (int j=0; j < M.ColNo(); j++)
      {
         printf("%.6f ", M(i,j));

      }
      printf("\n");
	}
}

//Creates an scaled square identity matrix
Matrix Make_eye(int NUM_RX, double Scale)
{
	int i,j;
	Matrix temp_eye(NUM_RX,NUM_RX);
	for(i=0;i<NUM_RX;i++)
	{
		for(j=0;j<NUM_RX;j++)
		{
			if(i==j)
				temp_eye(i,j) = 1 * Scale;
			else
				temp_eye(i,j) = 0;
		}
	}
	return temp_eye;
}

//Creates a channel matrix of dimension ACTUAL_NUM_RX by ACTUAL_NUM_TX
//the elements of the channel matrix are 0-mean variance-1/2 realizations
Matrix makeChannel(int ACTUAL_NUM_RX,int ACTUAL_NUM_TX)
{
		double x1 = 0;
		double x2 = 0;
		double y1,y2;
		int NUM_RX = 2 * ACTUAL_NUM_RX;
		int NUM_TX = 2 * ACTUAL_NUM_TX;
		Matrix H(NUM_RX,NUM_TX);
		double rt_2 = 1 / sqrt((double)2);
		int i,j;
		for (i=0; i < ACTUAL_NUM_RX; i++)
		{
			for (j=0; j < ACTUAL_NUM_TX; j++)
			{	
			
			    x1 = (double)rand()/RAND_MAX;
				x2 = (double)rand()/RAND_MAX;
				y1 = (double)sqrt( - 2 *log(x1) ) *cos( 2 *pi* x2 );
				y2 = (double)sqrt( - 2 *log(x1) ) *sin( 2 *pi* x2 );
				while( (x1 == 0) || (x2 == 0) || (y1 > 1000) || (y2 > 1000) 
					    || (y1 == 0) || (y2 == 0) || (y1 < -1000) || (y2 < -1000) )
				{
					x1 = (double)rand()/RAND_MAX;
					x2 = (double)rand()/RAND_MAX;
					y1 = (double)sqrt( - 2 *log(x1) ) *cos( 2 *pi* x2 );
					y2 = (double)sqrt( - 2 *log(x1) ) *sin( 2 *pi* x2 );
				}	
				H(i,j) = rt_2 * y1;
				H(i,j + ACTUAL_NUM_TX) = rt_2 * y2;
			}
		}
		for (i=0; i < ACTUAL_NUM_RX; i++)
		{
			for (j=0; j < ACTUAL_NUM_TX; j++)
			{	
				H(ACTUAL_NUM_RX + i, j) = -H(i, ACTUAL_NUM_TX + j);
				H(ACTUAL_NUM_RX + i, ACTUAL_NUM_TX + j) = H(i,  j);
			}
		}
		return H;
}

//Calculates the required scale factor for the transmitted symbols
// based on the required receive SNR and the number of transmit antennas
void getScale(double & scale, double SCALE, double & snr, double SNR,int ACTUAL_NUM_TX,int PAM)
{		
	if(SCALE == 0)
	{
		if(PAM == 2)
		  scale = 1/sqrt((double) 2);
		else if(PAM == 4)
		  scale = 1/sqrt((double) 10);
		else if(PAM == 8)
		  scale = 1/sqrt((double) 42);
		else if(PAM == 16)
		  scale = 1/sqrt((double) 170);
		else
			printf("ERROR num 2000");
	
		scale = scale*sqrt(SNR)/sqrt((double) ACTUAL_NUM_TX);
		snr = SNR;
	}
	else
	{
		if(PAM == 2)
		  snr = 1;
		else if(PAM == 4)
		  snr = 10;
		else if(PAM == 8)
		  snr = 42;
		else if(PAM == 16)
		  snr = 170;
		scale = SCALE;
		snr = snr * scale * scale * ACTUAL_NUM_TX;
	}
}

//Obtains a column vector of unscaled symbols in the PAM constellation
Matrix getInputs(int NUM_TX, int PAM)
{
		Matrix u(NUM_TX,1);
		double a;
		int b,i;

		for (i = 0; i < NUM_TX;i++)
		{
		  a = (double) 1000 * rand()/RAND_MAX ;
		  b = (int) a;
		  a = (double) 2*((b % PAM) - (((double)PAM - 1)/2));
		  u(i,0) = a;
		}
		return u;
}

//Creates a vector of noise realizations from a distribution
//  that has 0-mean variance-1/2
Matrix makeNoise(int NUM_RX)
{		
		int i;
		Matrix n(NUM_RX,1);
		int rx_odd = NUM_RX %2;
		double rt_2 = 1 / sqrt((double) 2);
		double x1,x2,y1,y2;
		for(i=0 ; i<NUM_RX/2 ; i++)
		{
			x1 = (double)rand()/RAND_MAX;
			x2 = (double)rand()/RAND_MAX;
			y1 = (double)sqrt( - 2 *log(x1) ) *cos( 2 *pi* x2 );
			y2 = (double)sqrt( - 2 *log(x1) ) *sin( 2 *pi* x2 );
			while( (x1 == 0) || (x2 == 0) || (y1 > 1000) || (y2 > 1000)
				    || (y1 == 0) || (y2 == 0) || (y1 < -1000) || (y2 < -1000) )
			{
				x1 = (double)rand()/RAND_MAX;
				x2 = (double)rand()/RAND_MAX;
				y1 = (double)sqrt( - 2 *log(x1) ) *cos( 2 *pi* x2 );
				y2 = (double)sqrt( - 2 *log(x1) ) *sin( 2 *pi* x2 );
			}	
			n(i,0) = rt_2 * y1;
			n(i+NUM_RX/2,0) = rt_2 * y2;
		}
		if(rx_odd)
		{
			x1 = (double)rand()/RAND_MAX;
			x2 = (double)rand()/RAND_MAX;
			y1 = (double)sqrt( - 2 *log(x1) ) *cos( 2 *pi* x2 );
			while( (x1 == 0) || (x2 == 0) || (y1 > 1000) || (y1 == 0) || (y1 < -1000) )
			{
				x1 = (double)rand()/RAND_MAX;
				x2 = (double)rand()/RAND_MAX;
				y1 = (double)sqrt( - 2 *log(x1) ) *cos( 2 *pi* x2 );
			}	
			n(NUM_RX-1,0) = rt_2 * y1;
		}
		return n;
}

//Counts how many bit errors a given symbol error equates to based on the 
//  grey coding of the bits and the PAM constellation
int count_bits(int PAM, double diff, double first, double second)
{
	int sum = 0;
	if(PAM == 2)
	{
		if(diff == 2 || diff == -2)
			sum = 1;
		else
			printf("WEIRD diff....\n\n");
	}
	else if(PAM == 4)
	{
		if(diff == 2 || diff == -2)
			sum = 1;
		else if(diff == 4 || diff == -4)
			sum = 2;
		else if(diff == 6 || diff == -6)
			sum = 1;
		else
			printf("WEIRD diff....PAM-4.. diff %f\n\n", diff);		
	}
	else if(PAM == 8)
	{
		if(diff == 2 || diff == -2)
			sum = 1;
		else if(diff == 4 || diff == -4)
			sum = 2;
		else if(diff == 6 || diff == -6)
		{
			if(first == -7 || first == -3 || first == 3 || first == 7 ||
				second == -7 || second == -3 || second == 3 || second == 7)
				sum = 1;
			else 
				sum = 3;
		}
		else if(diff == 8 || diff == -8)
			sum = 2;
		else if(diff == 10 || diff == -10)
		{
			if(first == -7 || first == -3 || first == 3 || first==7)
				sum = 3;
			else
				sum = 1;
		}
		else if(diff == 12 || diff == -12)
			sum = 2;
		else if(diff == 14 || diff == -14)
			sum = 1;
		else
			printf("WEIRD diff....\n\n");		
	}
	else if(PAM == 16)
	{
		if(diff == 2 || diff == -2)
			sum = 1;
		else if(diff == 4 || diff == -4)
			sum = 2;
		else if(diff == 6 || diff == -6)
		{
			double l_temp,u_temp;
			if(first < second)
			{
				l_temp = first;
				u_temp = second;
			}
			else
			{
				l_temp = second;
				u_temp = first;
			}
			if( (l_temp == -13 && u_temp == -7) || (l_temp == -9 && u_temp == -3)
				|| (l_temp == -5 && u_temp == 1) || (l_temp == -1 && u_temp == 5) 
				|| (l_temp == 3 && u_temp == 9) || (l_temp == 7 && u_temp == 13))
				sum = 3;
			else
				sum = 1;
		}
		else if(diff == 8 || diff == -8)
			sum = 2;
		else if(diff == 10 || diff == -10)
		{
			if(first == -13 || first == 13 || second == -13 || second == 13 || 
				(first == 5 && second == -5) || (first == -5 && second == 5))
				sum = 1;
			else
				sum = 3;
		}
		else if(diff == 12 || diff == -12)
		{
			if(first == -11 || second == -11 || first == -9 || second == -9 || 
				first == 9 || second == 9 || first == 11 || second == 11)
				sum = 4;
			else
				sum = 2;
		}
		else if(diff == 14 || diff == -14)
		{
			if(first == -7 || second == -7 || first == -15 || second == -15 || 
				first == 15 || second == 15)
				sum = 1;
			else
				sum = 3;
		}
		else if(diff == 16 || diff == -16)
			sum = 2;
		else if(diff == 18 || diff == -18)
		{
			if(first == 9 || second == 9)
				sum = 1;
			else
				sum = 3;
		}
		else if(diff == 20 || diff == -20)
		{
			if(first == 9 || second == 9 ||first == 11 || second == 11)
				sum = 2;
			else
				sum = 4;
		}
		else if(diff == 22 || diff == -22)
		{
			if(first == 11 || second == 11)
				sum = 1;
			else
				sum = 3;
		}
		else if(diff == 24 || diff == -24)
			sum = 2;
		else if(diff == 26 || diff == -26)
		{
			if(first == 13 || second == 13)
				sum = 1;
			else
				sum = 3;
		}
		else if(diff == 28 || diff == -28)
			sum = 2;
		else if(diff == 30 || diff == -30)
			sum = 1;
	}

	return sum;
}

//Finds the nearest constellation symbol for a vector
Matrix get_Eq_bits(int NUM_TX,Matrix Equalized,int PAM)
{
		Matrix Eq_bits(NUM_TX,1);
		int i;
		if(PAM == 2)
		{
		     for(i=0;i<NUM_TX;i++)
		     {
			if(Equalized(i,0) < 0)
				Eq_bits(i,0) = -1;
			else
				Eq_bits(i,0) = 1;
		     } 
		}
		else if(PAM == 4)
		{
		     for(i=0;i<NUM_TX;i++)
		     {
		        if(Equalized(i,0) < -2)
				Eq_bits(i,0) = -3;
			else if(Equalized(i,0) < 0)
				Eq_bits(i,0) = -1;
			else if(Equalized(i,0) < 2)
			  	Eq_bits(i,0) = 1;
			else
			        Eq_bits(i,0) = 3;
		     }

		}
		else if(PAM == 8)
		{
		     for(i=0;i<NUM_TX;i++)
		     {
			if(Equalized(i,0) < -6)
				Eq_bits(i,0) = -7;
			else if(Equalized(i,0) < -4)
				Eq_bits(i,0) = -5;
			else if(Equalized(i,0) < -2)
				Eq_bits(i,0) = -3;
			else if(Equalized(i,0) < 0)
				Eq_bits(i,0) = -1;
			else if(Equalized(i,0) < 2)
				Eq_bits(i,0) = 1;
			else if(Equalized(i,0) < 4)
				Eq_bits(i,0) = 3;
			else if(Equalized(i,0) < 6)
				Eq_bits(i,0) = 5;
			else
				Eq_bits(i,0) = 7;
		     }
		}
		else if(PAM == 16)
		{
		     for(i=0;i<NUM_TX;i++)
		     {			
		        if(Equalized(i,0) < -14)
				Eq_bits(i,0) = -15;
			else if(Equalized(i,0) < -12)
				Eq_bits(i,0) = -13;
			else if(Equalized(i,0) < -10)
				Eq_bits(i,0) = -11;
			else if(Equalized(i,0) < -8)
			        Eq_bits(i,0) = -9;
			else if(Equalized(i,0) < -6)
				Eq_bits(i,0) = -7;
			else if(Equalized(i,0) < -4)
				Eq_bits(i,0) = -5;
			else if(Equalized(i,0) < -2)
				Eq_bits(i,0) = -3;
			else if(Equalized(i,0) < 0)
				Eq_bits(i,0) = -1;
			else if(Equalized(i,0) < 2)
				Eq_bits(i,0) = 1;
			else if(Equalized(i,0) < 4)
				Eq_bits(i,0) = 3;
			else if(Equalized(i,0) < 6)
				Eq_bits(i,0) = 5;
			else if(Equalized(i,0) < 8)
				Eq_bits(i,0) = 7;			
			else if(Equalized(i,0) < 10)
				Eq_bits(i,0) = 9;
			else if(Equalized(i,0) < 12)
				Eq_bits(i,0) = 11;
			else if(Equalized(i,0) < 14)
				Eq_bits(i,0) = 13;
			else
				Eq_bits(i,0) = 15;
		     }
		}
		return Eq_bits;
}

//Used to find the row of W that has the smallest norm
int get_smallest(int NUM_TX, int NUM_RX, Matrix W, Matrix picked, int & up_to)
{
	double smallest_size = 99999999;
	int smallest_int = -1;
	int temp_up_to = -1;
	for(int i = 0; i < NUM_TX; i++)
	{
		if( picked(i,0) == 0)
		{
			temp_up_to++;
			Matrix temp(NUM_RX,1);
			for(int j = 0;j< NUM_RX; j++)
				temp(j,0) = W(temp_up_to,j);
			double size = eSize(temp,NUM_RX);
			if(size < smallest_size)
			{
				smallest_int = i;
				smallest_size = size;
				up_to = temp_up_to;
			}
		}
	}
	return smallest_int;
}			

//Implements the VBLAST decoding, a nulling and cancelling decoder
Matrix vblast_decoder(Matrix H, int NUM_TX, int NUM_RX, Matrix y, double scale, int PAM)
{
	Matrix y_VBLAST = y;
	Matrix picked(NUM_TX,1);
	Matrix VBLAST_OUT(NUM_TX,1);
	int i;
	for(i=0;i<NUM_TX;i++) picked(i,0) = 0;
	for(i=0;i<NUM_TX;i++)
	{
		int dim_left = NUM_TX - i;
		Matrix W(dim_left,NUM_RX);
		Matrix H_VBLAST(NUM_RX,dim_left);
		int upto = 0;
		for(int k = 0;k< NUM_TX; k++)
		{
			if(picked(k,0) == 0)
			{
				for(int j = 0;j < NUM_RX; j++)
					H_VBLAST(j,upto) = H(j,k);
				upto++;
			}
		}
		W = !( (~H_VBLAST) * H_VBLAST ) * (~H_VBLAST);
		int up_to = 0;
		int smallest_w = get_smallest(NUM_TX, NUM_RX, W, picked, up_to);
		picked(smallest_w,0) = 1;
		Matrix w_smallest(NUM_RX,1);
		Matrix h_smallest(NUM_RX,1);
		for(int j= 0;j< NUM_RX; j++)
		{
			w_smallest(j,0) = W(up_to,j);
			h_smallest(j,0) = H_VBLAST(j,up_to);
		}
		Matrix temp_VBLAST = (~w_smallest) * y_VBLAST;
		temp_VBLAST(0,0) = temp_VBLAST(0,0)/scale;
		Matrix VBLAST_sym = get_Eq_bits(1,temp_VBLAST,PAM);
		VBLAST_OUT(smallest_w,0) = VBLAST_sym(0,0);
		VBLAST_sym(0,0) = VBLAST_sym(0,0) * scale;
		y_VBLAST = y_VBLAST - ( h_smallest * VBLAST_sym );
	}
	return VBLAST_OUT;
}

Matrix makeV22_channel(Matrix H)
{
	Matrix tempH(8,8);

	tempH(0,0) = H(0,0);
	tempH(1,0) = 0;
	tempH(2,0) = H(0,2);
	tempH(3,0) = 0;
	tempH(4,0) = H(1,0);
	tempH(5,0) = 0;
	tempH(6,0) = H(1,2);
	tempH(7,0) = 0;

	tempH(0,1) = -H(0,2);
	tempH(1,1) = 0;
	tempH(2,1) = H(0,0);
	tempH(3,1) = 0;
	tempH(4,1) = -H(1,2);
	tempH(5,1) = 0;
	tempH(6,1) = H(1,0);
	tempH(7,1) = 0;

	tempH(0,2) = H(0,1);
	tempH(1,2) = 0;
	tempH(2,2) = H(0,3);
	tempH(3,2) = 0;
	tempH(4,2) = H(1,1);
	tempH(5,2) = 0;
	tempH(6,2) = H(1,3);
	tempH(7,2) = 0;

	tempH(0,3) = -H(0,3);
	tempH(1,3) = 0;
	tempH(2,3) = H(0,1);
	tempH(3,3) = 0;
	tempH(4,3) = -H(1,3);
	tempH(5,3) = 0;
	tempH(6,3) = H(1,1);
	tempH(7,3) = 0;

	tempH(0,4) = 0;
	tempH(1,4) = H(0,0);
	tempH(2,4) = 0;
	tempH(3,4) = H(0,2);
	tempH(4,4) = 0;
	tempH(5,4) = H(1,0);
	tempH(6,4) = 0;
	tempH(7,4) = H(1,2);

	tempH(0,5) = 0;
	tempH(1,5) = -H(0,2);
	tempH(2,5) = 0;
	tempH(3,5) = H(0,0);
	tempH(4,5) = 0;
	tempH(5,5) = -H(1,2);
	tempH(6,5) = 0;
	tempH(7,5) = H(1,0);

	tempH(0,6) = 0;
	tempH(1,6) = H(0,1);
	tempH(2,6) = 0;
	tempH(3,6) = H(0,3);
	tempH(4,6) = 0;
	tempH(5,6) = H(1,1);
	tempH(6,6) = 0;
	tempH(7,6) = H(1,3);

	tempH(0,7) = 0;
	tempH(1,7) = -H(0,3);
	tempH(2,7) = 0;
	tempH(3,7) = H(0,1);
	tempH(4,7) = 0;
	tempH(5,7) = -H(1,3);
	tempH(6,7) = 0;
	tempH(7,7) = H(1,1);

	return tempH;
}

Matrix makeLD22_channel(Matrix H)
{
	Matrix tempH(8,8);
	double rt = 1 / sqrt((double)2);

	tempH(0,0) = rt * H(0,0);
	tempH(1,0) = rt * H(0,1);
	tempH(2,0) = rt * H(0,2);
	tempH(3,0) = rt * H(0,3);
	tempH(4,0) = rt * H(1,0);
	tempH(5,0) = rt * H(1,1);
	tempH(6,0) = rt * H(1,2);
	tempH(7,0) = rt * H(1,3);

	tempH(0,1) = -rt * H(0,2);
	tempH(1,1) = -rt * H(0,3);
	tempH(2,1) = rt * H(0,0);
	tempH(3,1) = rt * H(0,1);
	tempH(4,1) = -rt * H(1,2);
	tempH(5,1) = -rt * H(1,3);
	tempH(6,1) = rt * H(1,0);
	tempH(7,1) = rt * H(1,1);

	tempH(0,2) = rt * H(0,1);
	tempH(1,2) = rt * H(0,0);
	tempH(2,2) = rt * H(0,3);
	tempH(3,2) = rt * H(0,2);
	tempH(4,2) = rt * H(1,1);
	tempH(5,2) = rt * H(1,0);
	tempH(6,2) = rt * H(1,3);
	tempH(7,2) = rt * H(1,2);

	tempH(0,3) = -rt * H(0,3);
	tempH(1,3) = -rt * H(0,2);
	tempH(2,3) = rt * H(0,1);
	tempH(3,3) = rt * H(0,0);
	tempH(4,3) = -rt * H(1,3);
	tempH(5,3) = -rt * H(1,2);
	tempH(6,3) = rt * H(1,1);
	tempH(7,3) = rt * H(1,0);

	tempH(0,4) = rt * H(0,0);
	tempH(1,4) = -rt * H(0,1);
	tempH(2,4) = rt * H(0,2);
	tempH(3,4) = -rt * H(0,3);
	tempH(4,4) = rt * H(1,0);
	tempH(5,4) = -rt * H(1,1);
	tempH(6,4) = rt * H(1,2);
	tempH(7,4) = -rt * H(1,3);

	tempH(0,5) = -rt * H(0,2);
	tempH(1,5) = rt * H(0,3);
	tempH(2,5) = rt * H(0,0);
	tempH(3,5) = -rt * H(0,1);
	tempH(4,5) = -rt * H(1,2);
	tempH(5,5) = rt * H(1,3);
	tempH(6,5) = rt * H(1,0);
	tempH(7,5) = -rt * H(1,1);

	tempH(0,6) = rt * H(0,1);
	tempH(1,6) = -rt * H(0,0);
	tempH(2,6) = rt * H(0,3);
	tempH(3,6) = -rt * H(0,2);
	tempH(4,6) = rt * H(1,1);
	tempH(5,6) = -rt * H(1,0);
	tempH(6,6) = rt * H(1,3);
	tempH(7,6) = -rt * H(1,2);

	tempH(0,7) = -rt * H(0,3);
	tempH(1,7) = rt * H(0,2);
	tempH(2,7) = rt * H(0,1);
	tempH(3,7) = -rt * H(0,0);
	tempH(4,7) = -rt * H(1,3);
	tempH(5,7) = rt * H(1,2);
	tempH(6,7) = rt * H(1,1);
	tempH(7,7) = -rt * H(1,0);


	return tempH;
}

Matrix addErrors( Matrix eH, int ACTUAL_NUM_RX, int ACTUAL_NUM_TX, double VAR_UNCERT)
{
	double x1 = 0;
	double x2 = 0;
	double y1,y2;
	int NUM_RX = 2 * ACTUAL_NUM_RX;
	int NUM_TX = 2 * ACTUAL_NUM_TX;
	Matrix H(NUM_RX,NUM_TX);
	double rt_2 = 1 / sqrt((double)2);
	double red = 0;
	if(VAR_UNCERT) red = sqrt(VAR_UNCERT);
	int i,j;
	for (i=0; i < ACTUAL_NUM_RX; i++)
	{
		for (j=0; j < ACTUAL_NUM_TX; j++)
		{	
		
			x1 = (double)rand()/RAND_MAX;
			x2 = (double)rand()/RAND_MAX;
			y1 = (double)sqrt( - 2 *log(x1) ) *cos( 2 *pi* x2 );
			y2 = (double)sqrt( - 2 *log(x1) ) *sin( 2 *pi* x2 );
			while( (x1 == 0) || (x2 == 0) || (y1 > 1000) || (y2 > 1000) 
					|| (y1 == 0) || (y2 == 0) || (y1 < -1000) || (y2 < -1000) )
			{
				x1 = (double)rand()/RAND_MAX;
				x2 = (double)rand()/RAND_MAX;
				y1 = (double)sqrt( - 2 *log(x1) ) *cos( 2 *pi* x2 );
				y2 = (double)sqrt( - 2 *log(x1) ) *sin( 2 *pi* x2 );
			}	
			H(i,j) = eH(i,j) + red * rt_2 * y1;
			H(i,j + ACTUAL_NUM_TX) = eH(i,j + ACTUAL_NUM_TX) + red * rt_2 * y2;
		}
	}
	for (i=0; i < ACTUAL_NUM_RX; i++)
	{
		for (j=0; j < ACTUAL_NUM_TX; j++)
		{	
			H(ACTUAL_NUM_RX + i, j) = -H(i, ACTUAL_NUM_TX + j);
			H(ACTUAL_NUM_RX + i, ACTUAL_NUM_TX + j) = H(i,  j);
		}
	}
	return H;
}
