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

#define CARRIERS 156
#define SYMBOLS 14
#define NUM_SUBFRAMES 1
#define WINDOW_SIZE 12

int mod_BPSK (int numbits, _Complex float *symbols);
int genRSsignalargerThan3RB(int u, int v, int m, int M_RS_SC, _Complex float *DMRSseq, int TxRxMode);
int largestprime_lower_than(int number);
int check_if_prime(int number);
void printZFCOEFF(_Complex float *correct, int nofcoeff);
void printCOEFF(int *correct, int nofcoeff);
int allocTest (int numbits, int *test);
void gridAllocationTest(int matrix[CARRIERS][SYMBOLS*NUM_SUBFRAMES], int *vector);
void gridAllocation(_Complex float grid[CARRIERS][SYMBOLS*NUM_SUBFRAMES], _Complex float *symbols);
void printMatrixTest(int matrix[CARRIERS][SYMBOLS*NUM_SUBFRAMES]);
void printMatrix(_Complex float matrix[CARRIERS][SYMBOLS*NUM_SUBFRAMES]);
void shiftWindowMatrix(int U[CARRIERS][WINDOW_SIZE]);
void shiftWindow(int *U);
void printMatrixWindow(int matrix[CARRIERS][WINDOW_SIZE]);
void putSamples(int U[CARRIERS][WINDOW_SIZE], int matrix[CARRIERS][SYMBOLS*NUM_SUBFRAMES], int column);
void averageVector(_Complex float *large, _Complex float *average);

//Modulation symbols and grid
_Complex float symbols[CARRIERS*SYMBOLS*NUM_SUBFRAMES];
_Complex float grid[CARRIERS][SYMBOLS*NUM_SUBFRAMES];

_Complex float average[SYMBOLS];

int test[CARRIERS*SYMBOLS*NUM_SUBFRAMES];
int testMatrix[CARRIERS][SYMBOLS*NUM_SUBFRAMES];

int U[CARRIERS][WINDOW_SIZE];

//DMRS Params
int M_RS_SC = 156;
int DMRS_length;
_Complex float DMRS_SEQ0[BUFFER_SZ];
_Complex float DMRS_SEQ1[BUFFER_SZ];

//Vectores de prueba
int W[CARRIERS][WINDOW_SIZE];

int Vector_Y[CARRIERS];
int error=0;
int deseada=5;
int eta=2;

int main() {
	clock_t t_ini, t_fin;
  	double secs;

  	t_ini = clock();

	//We modulate one subframe of samples
    /*mod_BPSK(CARRIERS*SYMBOLS*NUM_SUBFRAMES, symbols);
	printf("Vector: \n");
	printZFCOEFF(symbols, CARRIERS*SYMBOLS*NUM_SUBFRAMES);*/

	allocTest(CARRIERS*SYMBOLS*NUM_SUBFRAMES, test);
	//printf("Vector: \n");
    //printCOEFF(test, CARRIERS*SYMBOLS*NUM_SUBFRAMES);

	gridAllocationTest(testMatrix, test);
	//printf("Matrix: \n");
	//printMatrixTest(testMatrix);

	t_fin = clock();

  	secs = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;
  	printf("%.16g milisegundos\n", secs * 1000.0);

	/*averageVector(symbols, average);
	printf("Averaged vector: \n");
	printZFCOEFF(average, SYMBOLS);*/

	/*for(int i=0;i<14;i++){
		shiftWindow(U);
		U[0] = *(average+i);
		printCOEFF(U, WINDOW_SIZE);
	}*/
	
    int Y=0;
    int prueba=0;
    //Inicializamos vector de pesos
    for(int i=0; i<WINDOW_SIZE; i++){
		for(int j=0; j<CARRIERS; j++){
			W[i][j]=1+j;
		}
	}

	//printf("W: \n");
	//printMatrixWindow(W);

    for (int k = 0; k < 14; k++)
    {
        //Shift window
		shiftWindowMatrix(U);

		//Set new value
		putSamples(U, testMatrix, k);

		//printf("U: \n");
		//printMatrixWindow(U);

        //Realizamos la suma de coeficientes
		for(int i=0; i<CARRIERS; i++){
			for(int j=0; j<WINDOW_SIZE; j++){
				Vector_Y[i] += W[i][j]*U[i][j];
				//printf("Output %d: %d	", i, Vector_Y[i]);
			}
			//Calculamos error
			if(k==3||k==10){
				error=deseada-Vector_Y[i];
				//printf("Error: %d\n", error);

				//Actualizamos los pesos
				for(int j=0; j<WINDOW_SIZE; j++){
					//printf("U: \n");
					//printMatrixWindow(U);
					W[i][j] = W[i][j] + eta*error*U[i][j];
				}
			}
		}
		/*printf("\n");

		printf("\n");
		printf("Vector of outputs: \n");

		printf("Pesos: \n");
		printMatrixWindow(W);

    	printCOEFF(Vector_Y, CARRIERS);*/

        /*
        //Actualizamos los pesos
        for(int j=0; j<WINDOW_SIZE; j++)
			W[j] = W[j] + eta*error*U[j];

		printf("Pesos: \n");
        printCOEFF(W,12);

        Y=0;*/

		for(int i=0; i<CARRIERS; i++){
			Vector_Y[i] = 0;
		}
    }
   
    //printCOEFF(Vector_Y,12);

	/*printf("Window: \n");
	printMatrixWindow(U);*/

	//We allocate grid received into a matrix to equalize each subcarrier
	/*gridAllocation(grid, symbols);
	printMatrix(grid);*/

    DMRS_length = genRSsignalargerThan3RB(0, 1, 10, M_RS_SC, DMRS_SEQ0, 0);
    printf("DMRS 0 Length: %d\n", DMRS_length);
    DMRS_length = genRSsignalargerThan3RB(1, 1, 10, M_RS_SC, DMRS_SEQ1, 0);
    printf("DMRS 1 Length: %d\n", DMRS_length);

  	t_fin = clock();

  	secs = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;
  	printf("%.16g milisegundos\n", secs * 1000.0);
    return 0;
}

int allocTest (int numbits, int *test){
	int i,j;

	for(i=0; i<(SYMBOLS*NUM_SUBFRAMES); i++){
		for(j=0; j<CARRIERS; j++){
			*(test+CARRIERS*i+j)= i;
		}
	}
	return(i);
}

void printMatrixTest(int matrix[CARRIERS][SYMBOLS*NUM_SUBFRAMES]){
	int i, j; 
    for (i = 0; i < CARRIERS; i++){
		for (j = 0; j < (SYMBOLS*NUM_SUBFRAMES); j++) 
			printf("%d	", matrix[i][j]); 
		printf("\n");
	}
} 

void printMatrixWindow(int matrix[CARRIERS][WINDOW_SIZE]){
	int i, j; 
    for (i = 0; i < 3; i++){
		for (j = 0; j < WINDOW_SIZE; j++) 
			printf("%d	", matrix[i][j]); 
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

void printCOEFF(int *correct, int nofcoeff){
	int i;
	for(i=0; i<nofcoeff; i++){
		printf("%d,", *(correct+i));
	}
	printf("\n");
}

void gridAllocationTest(int matrix[CARRIERS][SYMBOLS*NUM_SUBFRAMES], int *vector){
	int i,j;
	for(i=0; i<(SYMBOLS*NUM_SUBFRAMES); i++){
		for(j=0; j<CARRIERS; j++){
			matrix[j][i] = vector[CARRIERS*i+j];
		}
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

void shiftWindowMatrix(int U[CARRIERS][WINDOW_SIZE]){
	int i,j;
	
	//Shift of vector
	for(i=WINDOW_SIZE-1; i > 0; i--){
		for(j=0; j<CARRIERS; j++){
			U[j][i] = U[j][i-1];
		}
	}
}

void shiftWindow(int *U){
	int i,j;
	//Shift of vector
	for(i=WINDOW_SIZE-1; i > 0; i--){
		U[i] = U[i-1];
	}
}

void putSamples(int U[CARRIERS][WINDOW_SIZE], int matrix[CARRIERS][SYMBOLS*NUM_SUBFRAMES], int column){
	int i,j;
	for(i=0; i<(SYMBOLS*NUM_SUBFRAMES); i++){
		for(j=0; j<CARRIERS; j++){
			U[j][0] = matrix[j][column];
		}
	}
}

void averageVector(_Complex float *large, _Complex float *average){
	int aux, avg;
	for(int i=0; i<SYMBOLS; i++){
		aux = 0;
		avg = 0;
		for(int j=0; j<CARRIERS; j++){
			aux += large[CARRIERS*i+j];
		}
		avg = aux/CARRIERS;
		average[i]=avg;
	}
}