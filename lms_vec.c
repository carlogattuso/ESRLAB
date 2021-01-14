#include <stdio.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define QAM4_LEVEL      0.7071    ///< QAM4 amplitude (RMS=1)
#define BUFFER_SZ	4096

#define CARRIERS 156
#define SYMBOLS 14
#define NUM_SUBFRAMES 1
#define WINDOW_SIZE 12

int mod_BPSK (int numbits, _Complex float *symbols);
void mod_TEST (int numbits, _Complex float *symbols);
void printZFCOEFF(_Complex float *correct, int nofcoeff);
void printVector(_Complex float *correct, int nofcoeff);
void shiftWindow(_Complex float *U);
int initWeights(_Complex float *W);

//Modulation symbols and grid
_Complex float symbols[BUFFER_SZ];
_Complex float equalized[BUFFER_SZ];
_Complex float U[BUFFER_SZ];
_Complex float W[BUFFER_SZ];
_Complex float Y = 0.0;
_Complex float error=0.0;
_Complex float desired=1.0;
float eta=0.01;

int main() {
	clock_t t_ini, t_fin;
  	double secs;

  	t_ini = clock();

    //We modulate one subframe of samples
    mod_BPSK(CARRIERS*SYMBOLS*NUM_SUBFRAMES, symbols);

    //Init vector of weights
    initWeights(W);

    for (int k = 0; k < SYMBOLS; k++){
        for(int i=0; i<CARRIERS; i++){
            
            //Shift window
		    shiftWindow(U);

            //Set new value to window
		    U[0] = symbols[CARRIERS*k+i];

            //We compute the output of the filter
            for(int j=0; j<WINDOW_SIZE; j++){
				Y += W[j]*U[j];
			}

            printf("Output: \n");
            if(__imag__ Y >= 0.0)printf("%f+%f*I,\n", __real__ Y,  __imag__ Y);
		    else printf("%f%f*I,\n", __real__ Y,  __imag__ Y);

            //Checking if it is a DMRS sequence
            //if(k==3||k==10) {
                //Calculate the error if it is a DMRS sample
                error = desired - Y;

                //We update weights accordingly
				for(int j=0; j<WINDOW_SIZE; j++) W[j] = W[j] + eta*error*U[j];
            //}

            //Append output to equalized vector
            equalized[CARRIERS*k+i] = Y;

            //Init output
            Y=0.0;
		}
    }

    printf("Equalized symbols: \n");
	printZFCOEFF(equalized, CARRIERS*SYMBOLS*NUM_SUBFRAMES);

    t_fin = clock();

  	secs = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;
  	printf("%.16g milisegundos\n", secs * 1000.0);
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