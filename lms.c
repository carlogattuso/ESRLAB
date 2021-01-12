#include <stdio.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <time.h>

#define QAM4_LEVEL      0.7071    ///< QAM4 amplitude (RMS=1)
#define BUFFER_SZ	2048
#define MOD_NO_DMRS_LENGTH 432
#define MOD_DMRS_LENGTH 504

#define CARRIERS 10
#define SYMBOLS 14
#define NUM_SUBFRAMES 1
#define WINDOW_SIZE 7

int mod_BPSK (int numbits, _Complex float *symbols);
int genRSsignalargerThan3RB(int u, int v, int m, int M_RS_SC, _Complex float *DMRSseq, int TxRxMode);
int largestprime_lower_than(int number);
int check_if_prime(int number);
void printZFCOEFF(_Complex float *correct, int nofcoeff);
void gridAllocation(_Complex float grid[CARRIERS][SYMBOLS*NUM_SUBFRAMES], _Complex float *symbols);
void printMatrix(_Complex float matrix[CARRIERS][SYMBOLS*NUM_SUBFRAMES]);
void printMatrixWindowComplex(_Complex float matrix[CARRIERS][WINDOW_SIZE]);
void shiftWindowMatrixComplex(_Complex float U[CARRIERS][WINDOW_SIZE]);
void putSamplesComplex(_Complex float U[CARRIERS][WINDOW_SIZE], _Complex float matrix[CARRIERS][SYMBOLS*NUM_SUBFRAMES], int column);
int initWeights(_Complex float *Weights);

//Modulation symbols and grid
_Complex float symbols[CARRIERS*SYMBOLS*NUM_SUBFRAMES];
_Complex float equalized[CARRIERS*SYMBOLS*NUM_SUBFRAMES];
_Complex float grid[CARRIERS][SYMBOLS*NUM_SUBFRAMES];

//Equalizer params
_Complex float U[CARRIERS][WINDOW_SIZE];
_Complex float weights[WINDOW_SIZE];
_Complex float W[CARRIERS][WINDOW_SIZE];
_Complex float Y = 0.0;
_Complex float error = 0.0;
_Complex float deseada = 5.0;
float eta = 0.005;

//DMRS Params
int M_RS_SC = 156;
int DMRS_length;
_Complex float DMRS_SEQ0[BUFFER_SZ];
_Complex float DMRS_SEQ1[BUFFER_SZ];
_Complex float DMRS_TX_TIME0[BUFFER_SZ];
_Complex float DMRS_RX_TIME0[BUFFER_SZ];
_Complex float DMRS_TX_TIME1[BUFFER_SZ];
_Complex float DMRS_RX_TIME1[BUFFER_SZ];

int main() {
	clock_t t_ini, t_fin;
  	double secs;

  	t_ini = clock();

	//We modulate one subframe of samples
    mod_BPSK(CARRIERS*SYMBOLS*NUM_SUBFRAMES, symbols);
	//printf("Vector: \n");
	//printZFCOEFF(symbols, CARRIERS*SYMBOLS*NUM_SUBFRAMES);

	initWeights(weights);
	printf("Init weights: \n");
	printZFCOEFF(weights, WINDOW_SIZE);

	//We allocate grid received into a matrix to equalize each subcarrier
	gridAllocation(grid, symbols);
	//printMatrix(grid);
	
    //Inicializamos vector de pesos
    for(int i=0; i<CARRIERS; i++){
		for(int j=0; j<WINDOW_SIZE; j++){
			W[i][j]=weights[j];
		}
	}

	printf("W: \n");
	printMatrixWindowComplex(W);

    for (int k = 0; k < SYMBOLS; k++)
    {
        //Shift window
		shiftWindowMatrixComplex(U);

		//Set new value
		putSamplesComplex(U, grid, k);

		//printf("U: \n");
		//printMatrixWindowComplex(U);

        //Realizamos la suma de coeficientes
		for(int i=0; i<CARRIERS; i++){
			Y=0.0;
			for(int j=0; j<WINDOW_SIZE; j++){
				Y += W[i][j]*U[i][j];
			}
			
			//printf("\n");
			//printf("Y: %f+%f*I \n", __real__ Y,  __imag__ Y);
			//printf("\n");

			//Añadimos la muestra al vector
			equalized[k*CARRIERS+i] = Y;

			//Calculamos error
			if(k==3||k==10){
				error=deseada-Y;
				//printf("Error: %f+%f*I \n", __real__ error,  __imag__ error);
				//printf("\n");
				
				//Actualizamos los pesos
				for(int j=0; j<WINDOW_SIZE; j++){
					//printf("U: %f+%f*I \n", __real__ U[i][j],  __imag__ U[i][j]);
					W[i][j] = W[i][j] + eta*error*U[i][j];
				}
			}
		}
		//printf("\n");
		//printf("Pesos: \n");
		//printMatrixWindowComplex(W);

        Y=0.0;
    }
   
	//printf("Vector: \n");
	//printZFCOEFF(equalized, CARRIERS*SYMBOLS*NUM_SUBFRAMES);

    DMRS_length = genRSsignalargerThan3RB(0, 1, 10, M_RS_SC, DMRS_SEQ0, 0);
    printf("DMRS 0 Length: %d\n", DMRS_length);
	//printZFCOEFF(DMRS_SEQ0, 156);
    DMRS_length = genRSsignalargerThan3RB(1, 1, 10, M_RS_SC, DMRS_SEQ1, 0);
    printf("DMRS 1 Length: %d\n", DMRS_length);

  	t_fin = clock();

  	secs = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;
  	printf("%.16g milisegundos\n", secs * 1000.0);
    return 0;
}

void printMatrixWindowComplex(_Complex float matrix[CARRIERS][WINDOW_SIZE]){
	int i, j; 
    for (i = 0; i < CARRIERS; i++){
		for (j = 0; j < WINDOW_SIZE; j++) 
			if(__imag__ matrix[i][j] >= 0.0)printf("%f+%f*I,	", __real__ matrix[i][j],  __imag__ matrix[i][j]);
			else printf("%f%f*I,	", __real__ matrix[i][j],  __imag__ matrix[i][j]);
		printf("\n");
	}
} 

void printMatrix(_Complex float matrix[CARRIERS][SYMBOLS*NUM_SUBFRAMES]){
	int i, j; 
    for (i = 0; i < CARRIERS; i++){
		for (j = 0; j < (SYMBOLS*NUM_SUBFRAMES); j++) 
			if(__imag__ matrix[i][j] >= 0.0)printf("%f+%f*I,	", __real__ matrix[i][j],  __imag__ matrix[i][j]);
			else printf("%f%f*I,	", __real__ matrix[i][j],  __imag__ matrix[i][j]);
		printf("\n");
	}
} 

int mod_BPSK (int numbits, _Complex float *symbols){
	int i, bit;

    bit = 0;
	for (i=0;i<numbits;i++){
		if(bit == 0){
            *(symbols+i)= QAM4_LEVEL + QAM4_LEVEL*I;
            bit = 1;
        }
		else {
            *(symbols+i)= -(QAM4_LEVEL + QAM4_LEVEL*I);
            bit = 0;
        }
	}
	return(i);
}

int genRSsignalargerThan3RB(int u, int v, int m, int M_RS_SC, _Complex float *DMRSseq, int TxRxMode){
	int i, sequence_length;
	double arg, rootidx;
	int N_RS_ZC = largestprime_lower_than(M_RS_SC);
	float qNfloat = ((float)(N_RS_ZC*(u+1)))/31.0;
	int qNint = qNfloat;
	float q = (float)(int)((qNfloat + 0.5)) + (float)v*(float)pow((double)(-1.0), (double)(int)(2.0*qNfloat));

	sequence_length=M_RS_SC;

	for(i=0; i<sequence_length; i++){
		arg=(((double)q)*((double)(i%N_RS_ZC)*((double)(i%N_RS_ZC)+1.0)))/(double)N_RS_ZC;
		__real__ DMRSseq[i]=(float)cos(arg)*5.0;
		__imag__ DMRSseq[i]=(float)sin(arg)*5.0;
	}
	return sequence_length;
}

int largestprime_lower_than(int number){

	int i; 

	for(i=1; i<number; i++){
		if(check_if_prime(number-i)==1) break;
	}
	return(number-i);
}

int check_if_prime(int number){

    int n, i, flag = 1;

    for(i = 2; i <= number/2; ++i)
    {
        // condition for nonprime number
        if(number%i == 0)
        {
            flag = 0;
            break;
        }
    }
	// flag 1 if prime
	// flag 0 if not
	return flag;
	
}

void printZFCOEFF(_Complex float *correct, int nofcoeff){
	int i;
	for(i=0; i<nofcoeff; i++){
		if(__imag__ *(correct+i) >= 0.0)printf("%f+%f*I,\n", __real__ *(correct+i),  __imag__ *(correct+i));
		else printf("%f%f*I,\n", __real__ *(correct+i),  __imag__ *(correct+i));
	}
}

void gridAllocation(_Complex float grid[CARRIERS][SYMBOLS*NUM_SUBFRAMES], _Complex float *symbols){
	int i,j;
	for(i=0; i<(SYMBOLS*NUM_SUBFRAMES); i++){
		for(j=0; j<CARRIERS; j++){
			grid[j][i] = symbols[CARRIERS*i+j];
		}
	}
}

void shiftWindowMatrixComplex(_Complex float U[CARRIERS][WINDOW_SIZE]){
	int i,j;
	
	//Shift of vector
	for(i=WINDOW_SIZE-1; i > 0; i--){
		for(j=0; j<CARRIERS; j++){
			U[j][i] = U[j][i-1];
		}
	}
}

void putSamplesComplex(_Complex float U[CARRIERS][WINDOW_SIZE], _Complex float matrix[CARRIERS][SYMBOLS*NUM_SUBFRAMES], int column){
	int i,j;
	for(i=0; i<(SYMBOLS*NUM_SUBFRAMES); i++){
		for(j=0; j<CARRIERS; j++){
			U[j][0] = matrix[j][column];
		}
	}
}

int initWeights(_Complex float *Weights){
	//We choose normally distributed init weights to begin equalizing
	*(Weights)=0.6604 + 0.0000*I;
	*(Weights+1)=-0.4210 + 0.0000*I;
  	*(Weights+2)=-0.7330 + 0.0000*I;
  	*(Weights+3)=0.3788 - 0.0000*I;
  	*(Weights+4)=0.2327 - 0.0000*I;
  	*(Weights+5)=-1.7549 - 0.0000*I;
  	*(Weights+6)=-1.4092 - 0.0000*I;
  	
	return(7);
}