void get_list(int NUM_TX,int PAM,int ui,int Li,double Si,
			  Matrix & Sorted,int & number, int dim)
{
	
	number = Li - ui;
	number = (int) number / 2;
	number = number + 1;
	double smallest = 999;
	int smallest_i = 0;
	int i = 0;
	double diff = 0;
	for(i = ui; i <= Li; i = i + 2)
	{
		diff = (double) i - Si;
		if( abs(diff) < abs(smallest) )
		{
			smallest = diff;
			smallest_i = i;
		}
	}
	int step = 2;
	if(smallest < 0)
		step = -2;
	int current_i = smallest_i;
	int flag = 0;
	int last_i = 0;
	int upto = 0;
	for(i = 0; i < number; i++)
	{
		last_i = current_i;
		upto = i;
		current_i = current_i + i * step;
		if(current_i > Li)
			flag = 1;
		if(current_i < ui)
			flag = 2;
		if(flag) break;
		Sorted(dim, i ) = (double) current_i;  
		step = -step;
	}
	if(flag)
	{
		current_i = last_i;
		int added = 0;
		if(flag == 1) added = -2; else added = 2;
		for(i = upto; i< number; i++)
		{
			current_i = current_i + added;
			Sorted(dim, i ) = (double) current_i; 
		}
	}
} 

Matrix se_sphere_decoder(Matrix H, int NUM_TX, int NUM_RX, Matrix Equalized, 
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
	int N[50];
	int x[50];

	double ep[50];
	int best_point[50];
	for(j= 0; j< NUM_TX; j++) best_point[j] = 0;
	L[i] = floor( sqrt( T(i,0)/Q(i,i) ) + S(i,0) );
	u[i] = ceil( -sqrt( T(i,0)/Q(i,i) ) + S(i,0) );
	if(L[i]%2 == 0) L[i]--;
	if(u[i]%2 == 0) u[i]++;
	adj_uL(PAM,u[i],L[i]);
	Matrix Sorted(NUM_TX,2*PAM);
	for(int ll = 0; ll < NUM_TX; ll++)
	{
		for(int pp = 0; pp < 2*PAM; pp++)
			Sorted(ll,pp) = 0;
	}
	get_list(NUM_TX,PAM,u[i],L[i],S(i,0),Sorted,N[i],i);
	x[i] = -1;
	tester = 0;
	while (1)
	{
		tester++;
		x[i]++;
		if(x[i] >= N[i])  //case: no more symbols at the layer
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
			ep[i] = Equalized(i,0) - Sorted(i,x[i]);
			T(i-1,0) = T(i,0) - Q(i,i) * ( S(i,0) - Sorted(i,x[i]) ) * 
				       ( S(i,0) - Sorted(i,x[i]) );
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
			get_list(NUM_TX,PAM,u[i],L[i],S(i,0),Sorted,N[i],i);
			x[i] = -1;
		}
		else   //case: lowest layer, see if the symbol has the smallest distance
		{
			double d_hat2 = T(NUM_TX-1,0) - T(0,0) + Q(0,0) * 
				            ( S(0,0) - Sorted(0,x[0]) ) * ( S(0,0) - Sorted(0,x[0]) );
			if(d_hat2 < d_2)  //case: smallest distance so far  
			{
				d_2 = d_hat2;
				T(NUM_TX-1,0) = d_2;
				int k;
				if(x[NUM_TX-1] < 0) continue;
				for(k = NUM_TX - 1; k > 0 ; k--)
					T(k-1,0) = T(k,0) - Q(k,k) * (S(k,0) - Sorted(k,x[k]) ) * 
					           (S(k,0) - Sorted(k,x[k]) );
				for(k = 0 ; k < NUM_TX; k++)
				{
					best_point[k] = Sorted(k,x[k]);
					L[k] = floor( sqrt( T(k,0)/Q(k,k) ) + S(k,0) );
					u[k] = ceil( -sqrt( T(k,0)/Q(k,k) ) + S(k,0) );
					if(L[k]%2 == 0) L[k]--;
					if(u[k]%2 == 0) u[k]++;
					adj_uL(PAM,u[k],L[k]);
					get_list(NUM_TX,PAM,u[k],L[k],S(k,0),Sorted,N[k],k);
				}
			}
		}
	} //END OF while(1)
	Matrix sphereBest(NUM_TX,1);
	for(j= 0; j< NUM_TX; j++) 
		sphereBest(j,0) = best_point[j];
	return sphereBest;
}

