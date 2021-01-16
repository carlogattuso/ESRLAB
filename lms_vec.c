#include <stdio.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define QAM4_LEVEL      0.7071    ///< QAM4 amplitude (RMS=1)
#define BUFFER_SZ	2048

#define CARRIERS 5
#define SYMBOLS 14
#define NUM_SUBFRAMES 1
#define WINDOW_SIZE 12

int mod_BPSK (int numbits, _Complex float *symbols);
void mod_TEST (int numbits, _Complex float *symbols);
void printZFCOEFF(_Complex float *correct, int nofcoeff);
void printVector(_Complex float *correct, int nofcoeff);
void shiftWindow(_Complex float *U);
int initWeights(_Complex float *W);
int addNoise (int numbits, _Complex float *symbols, _Complex float *symbols_noise);

//Modulation symbols and grid
_Complex float symbols[BUFFER_SZ];
_Complex float symbols_noise[BUFFER_SZ];
_Complex float equalized[BUFFER_SZ];
_Complex float U[BUFFER_SZ];
_Complex float W[BUFFER_SZ];
_Complex float Y = 0.0;
_Complex float error=0.0;
_Complex float desired=1.0;
float eta=0.05;

int main() {
	clock_t t_ini, t_fin;
  	double secs;

	srand(time(0));
  	t_ini = clock();

    //We modulate one subframe of samples
    mod_BPSK(CARRIERS*SYMBOLS*NUM_SUBFRAMES, symbols);
	addNoise(CARRIERS*SYMBOLS*NUM_SUBFRAMES, symbols, symbols_noise);

	float test;

	for(int i=0;i<20;i++){
		test = (float) (rand() % 6) - 3.0;
	}

	printf("Test number: %f\n", test);

    //Init vector of weights
    initWeights(W);

    for (int k = 0; k < SYMBOLS; k++){
        for(int i=0; i<CARRIERS; i++){
            
			//Init output
            Y=0.0;

            //Shift window
		    shiftWindow(U);

            //Set new value to window
			U[0] = symbols[CARRIERS*k+i];
		    //U[0] = symbols_noise[CARRIERS*k+i];

            //We compute the output of the filter
            for(int j=0; j<WINDOW_SIZE; j++){
				Y += W[j]*U[j];
			}

            //Checking if it is a DMRS sequence
            //if(k==3||k==10) {
                //Calculate the error if it is a DMRS sample
				desired = symbols[CARRIERS*k+i];
                error = desired - Y;

				printf("\n");
				printf("CONJUGATED: %f+%f*I \n", __real__ error,  __imag__ error);
				printf("\n");
				
                //We update weights accordingly
				for(int j=0; j<WINDOW_SIZE; j++) W[j] = W[j] + eta*conj(error)*U[j];
            //}

            //Append output to equalized vector
            equalized[CARRIERS*k+i] = Y;
		}
    }

    t_fin = clock();

  	secs = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;
  	printf("%.16g milisegundos\n", secs * 1000.0);
}

int mod_BPSK (int numbits, _Complex float *symbols){
	int i, bit=0, sign =0;

	for (i=0;i<numbits;i++){
		bit = rand() % 2;
		if(bit == 1){
            /**(symbols+i)= (QAM4_LEVEL + QAM4_LEVEL*I);*/
			*(symbols+i)= (1.0 + 1.0*I);
        }
		else {
            /**(symbols+i)= -(QAM4_LEVEL + QAM4_LEVEL*I);*/
			*(symbols+i)= -(1.0 + 1.0*I);
        }
	}
	return(i);
}

int addNoise (int numbits, _Complex float *symbols, _Complex float *symbols_noise){
	int i, bit;
	float test;

	for (i=0;i<numbits;i++){
		test = (float) ((rand() % 6) - 3.0)/10;
		symbols_noise[i] = symbols[i]+test*I;
	}
	return(i);
}

void mod_TEST (int numbits, _Complex float *symbols){
	int i, j;

    for (int k = 0; k < SYMBOLS; k++)
        for(int i=0; i<CARRIERS; i++)
            *(symbols+CARRIERS*k+i)= (float) (k+1) + (float) (k+1)*I;
}


void printZFCOEFF(_Complex float *correct, int nofcoeff){
	int i;
	for(i=0; i<nofcoeff; i++){
		if(__imag__ *(correct+i) >= 0.0)printf("%f+%f*I,\n", __real__ *(correct+i),  __imag__ *(correct+i));
		else printf("%f%f*I,\n", __real__ *(correct+i),  __imag__ *(correct+i));
	}
}

void printVector(_Complex float *correct, int nofcoeff){
	int i;
    
	for(i=0; i<nofcoeff; i++){
		if(__imag__ *(correct+i) >= 0.0) printf("%f+%f*I,    ", __real__ *(correct+i),  __imag__ *(correct+i));
		else printf("%f%f*I,    ", __real__ *(correct+i),  __imag__ *(correct+i));
	}
    printf("\n");
}

void shiftWindow(_Complex float *U){
	for(int i=WINDOW_SIZE-1; i > 0; i--) U[i] = U[i-1];
}

int initWeights(_Complex float *Weights){
	//We choose normally distributed init weights to begin equalizing
	*(Weights)=0.0000 + 0.0000*I;
	*(Weights+1)=0.0000 + 0.0000*I;
  	*(Weights+2)=0.0000 + 0.0000*I;
  	*(Weights+3)=0.0000 + 0.0000*I;
  	*(Weights+4)=0.0000 + 0.0000*I;
  	*(Weights+5)=0.0000 + 0.0000*I;
  	*(Weights+6)=0.0000 + 0.0000*I;
  	
	return(7);
}