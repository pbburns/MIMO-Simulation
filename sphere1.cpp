/////////////////////////////////////////////////////
// file:        sphere1.h                          //
// copywrite:   Philippe Bergeron-Burns            //
// description: Implements the regular SD          //
/////////////////////////////////////////////////////

void adj_uL( int PAM, int & u, int & L);

//////////////////////////////////////////////////////
//Inputs  H          - the channel matrix
//        NUM_TX     - the number of transmit antennas
//        NUM_RX     - the number of receive antennas
//        Equalized  - the received equalized vector
//        difference - the difference between Equalized and Eq_bits
//        PAM        - the symbol constellation
//Outputs tester     - # times the SD goes through the main loop
//        return (best_bits) - the vector of symbols that are ML
//////////////////////////////////////////////////////
Matrix sphere_decoder(Matrix H, int NUM_TX, int NUM_RX, Matrix Equalized, 
					  Matrix difference, int PAM, int & tester)
{
	Matrix QQ(NUM_RX,NUM_TX);
	Matrix R(NUM_TX,NUM_TX);
	
	//now do the QR factorization of H
	QR(QQ,R,H,NUM_TX,NUM_RX);

	Matrix Q(NUM_TX, NUM_TX);
	int j,k;
	for(j = 0; j< NUM_TX; j++)
	{
		for(k = 0; k< NUM_TX; k++)
		{
			if(k < j) Q(j,k) = 0;
			else if(j == k) Q(j,j) = R(j,j) * R(j,j);
			else Q(j,k) = R(j,k) / R(j,j);
		}
	}

	Matrix prod = R * difference;
	double d_2 = eSize(prod, NUM_TX);
	d_2 = d_2 * d_2 * 1.01;
	Matrix S = Equalized;
	Matrix T(NUM_TX,1);
	for(j = 0; j < NUM_TX; j++) T(j,0) = 0;

	T(NUM_TX-1,0) = d_2;
	int i = NUM_TX-1;

	//won't work if NUM_TX is more than 50
	int L[50];
	int u[50];

	double ep[50];
	int best_point[50];
	for(j= 0; j< NUM_TX; j++) best_point[j] = 0;
	L[i] = floor( sqrt( T(i,0)/Q(i,i) ) + S(i,0) );
	u[i] = ceil( -sqrt( T(i,0)/Q(i,i) ) + S(i,0) );
	if(L[i]%2 == 0) L[i]--;
	if(u[i]%2 == 0) u[i]++;
	adj_uL(PAM,u[i],L[i]);
	u[i] = u[i] - 2;
	tester = 0;
	while (1)
	{
		tester++;
		u[i] = u[i] + 2;
		if(u[i] > L[i])  //case: no more symbols at the layer
		{
			if(i == NUM_TX-1)
			{
				break;
			}
			i = i + 1;
			continue;
		}
		if(i>0)  //normal case for a symbol: move down to a lower layer
		{
			ep[i] = Equalized(i,0) - u[i];
			T(i-1,0) = T(i,0) - Q(i,i) * ( S(i,0) - u[i] ) * ( S(i,0) - u[i] );
			double sum_t = 0;
			for(j = i ; j < NUM_TX; j++)
				sum_t += Q(i-1,j) * ep[j];
			S(i-1,0) = Equalized(i-1,0) + sum_t;
			i = i - 1;
			L[i] = floor( sqrt( T(i,0)/Q(i,i) ) + S(i,0) );
			u[i] = ceil( -sqrt( T(i,0)/Q(i,i) ) + S(i,0) );
			if(L[i]%2 == 0) L[i]--;
			if(u[i]%2 == 0) u[i]++;
			adj_uL(PAM,u[i],L[i]);
			u[i] = u[i] - 2;
		}
		else   //case: lowest layer, see if the symbol has the smallest distance
		{
			double d_hat2 = T(NUM_TX-1,0) - T(0,0) + Q(0,0) 
				              * ( S(0,0) - u[0] ) * ( S(0,0) - u[0] );
			if(d_hat2 < d_2)  //case: smallest distance so far  
			{
				d_2 = d_hat2;
				for(j= 0; j< NUM_TX; j++) 
					best_point[j] = u[j];
				T(NUM_TX-1,0) = d_2;
				i = NUM_TX-1;
				L[i] = floor( sqrt( T(i,0)/Q(i,i) ) + S(i,0) );
				u[i] = ceil( -sqrt( T(i,0)/Q(i,i) ) + S(i,0) );
				if(L[i]%2 == 0) L[i]--;
				if(u[i]%2 == 0) u[i]++;
				adj_uL(PAM,u[i],L[i]);
				u[i] = u[i] - 2;
			}
		}
	} //END OF while(1)
	Matrix sphereBest(NUM_TX,1);
	for(j= 0; j< NUM_TX; j++) 
		sphereBest(j,0) = best_point[j];
	return sphereBest;
}

void adj_uL( int PAM, int & u, int & L)
{
	if(PAM)
	{
		int t = PAM-1;
		if(u < -t) u = -t;
		if(L > t) L = t;
	}
}
