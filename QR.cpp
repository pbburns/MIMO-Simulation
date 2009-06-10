/////////////////////////////////////////////////////
// file:        QR.h                               //
// copywrite:   Philippe Bergeron-Burns            //
// description: Does QR factorization of an        //
//              arbitrary matrix                   //
/////////////////////////////////////////////////////

void updateQR(Matrix & Q, Matrix & R, int NUM_TX, int NUM_RX,int i);

//////////////////////////////////////////////////////
//Inputs  H          - the channel matrix
//        NUM_TX     - the number of transmit antennas
//        NUM_RX     - the number of receive antennas
//Outputs Q          - The unique matrix Q such that QtQ = I and H = QR
//        R          - The upper triangular matrix where H = QR
///////////////////////////////////////////////////////
void QR(Matrix & Q, Matrix & R, Matrix H, int NUM_TX, int NUM_RX)
{
	Matrix fullR(NUM_RX,NUM_TX);
	Matrix fullQ(NUM_RX,NUM_RX);	
	
	//initialize the fullQ matrix to I and the fullR matrix to H
	int i,j; 
	for(i = 0; i < NUM_RX; i++)
	{
		for(j = 0 ; j < NUM_RX; j++)
		{
			if(i == j) fullQ(i,j) = 1;
			else fullQ(i,j) = 0;
		}
		for(j = 0; j< NUM_TX; j++)
			fullR(i,j) = H(i,j);		
	}

	//now using the householder method do a QR factorization
	for(i = 0; i < NUM_TX; i++)
	{
		updateQR(fullQ,fullR,NUM_TX,NUM_RX,i);
	}

	//convert the fullR matrix back to R and the
	//            fullQ matrix back to Q
	for(i = 0; i < NUM_RX; i++)
		for(j = 0; j< NUM_TX; j++)
			Q(i,j) = fullQ(i,j);		
	for(i = 0; i < NUM_TX; i++)
		for(j = 0; j< NUM_TX; j++)
			R(i,j) = fullR(i,j);		

	
}

void updateQR(Matrix & Q, Matrix & R, int NUM_TX, int NUM_RX,int i)
{
		
		int ni = NUM_RX - i;
		int ti = NUM_TX - i;
		Matrix Rtemp(ni,ti);
		Matrix Qtemp(ni,ni);
		Matrix Atemp(ni,ti);
		Matrix Itemp(ni,ni);
		int k,l;
		for(k = i; k < NUM_RX; k++)
		{
			for(l = i ; l < NUM_TX; l++)
			{
				Atemp(k-i , l-i) = R(k,l);
			}
		}
		for(k = 0; k < ni; k++)
		{
			for(l = 0; l < ni; l++)
			{
				if(k == l) Itemp(k,l) = 1;
				else Itemp(k,l) = 0;
			}
		}
		Matrix x(ni,1);
		for(k = 0; k < ni; k++)
			x(k,0) = Atemp(k,0);
		double xSize = eSize(x,ni);
		int sign; 
		if ( x(0,0) >= 0 ) sign = 1; 
		else sign = -1;
		Matrix g(ni,1);
		g(0,0) = x(0,0) + sign * xSize;
		for(k = 1; k < ni; k++)
			g(k,0) = x(k,0);
		double gSize = eSize(g,ni);
		double Qscale = 2/(gSize * gSize);
		Matrix ggT(ni,ni);
		ggT = g * (~g);
		for(k = 0; k < ni; k++)
			for(l = 0; l < ni; l++)
				ggT(k,l) = Qscale * ggT(k,l);
		Qtemp = Itemp - ggT;
		Rtemp = Qtemp * Atemp;

		Matrix Qeff(NUM_RX,NUM_RX);
		for(k = 0; k < NUM_RX; k++)
			for(l = 0; l < NUM_RX; l++)
				Qeff(k,l) = 0;
		for(k = 0; k < i; k++) Qeff(k,k) = 1;
		
		for(k = i; k < NUM_RX; k++)
		{
			for(l = i; l < NUM_TX; l++)
				R(k,l) = Rtemp(k-i,l-i);
			for(l = i; l < NUM_RX; l++)
				Qeff(k,l) = Qtemp(k-i,l-i);
		}
		Q = Q * Qeff; 
}

