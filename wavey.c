//config
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define FPS 20 //Frame rate of visual, also my sample rate of audio data
#define s_factor 10 //Factor by which to smoothen the data. (lose hi end freqs)


// ### Fourier Transform ###


#define DFT -1  //Sign in exponent; determines if transform is forward or inverse.
#define IDFT 1

//complex number struct
typedef struct {
	float re;
	float im; 
} cnum;

//cnum magnitude function
float mag(cnum z) {
	return sqrt(z.re * z.re + z.im * z.im);
}

//multiply two cnums
cnum multiply(cnum w, cnum z) {
	return (cnum){w.re*z.re - w.im*z.im,w.re*z.im + w.im*z.re};
}

//e^{+/- 2*pi*n*k/N}
cnum factor(int n, int k, int N, int sign) {
	return (cnum){cos(2*M_PI*n*k/N),sin(sign*2*M_PI*n*k/N)};
}


//DFT or IDFT depending on the first parameter
cnum *transform(int sign, cnum *hs, int N) {

	cnum *results = (cnum *)malloc((size_t)N * sizeof(cnum));
	if (results == NULL) {
		printf("\nProblem with memory allocation\n\n");
		exit(EXIT_FAILURE);
	}

	//Iterate for each term in the transform
	for (int n=0;n<N;n++) {
		cnum h = {0,0};
		
		//Iterate for the sum in the transform
		for (int k=0;k<N;k++) {
			cnum c = multiply(hs[k],factor(n,k,N,sign));

			//If IDFT, normalise results.
			if (sign == 1) {
				h.re += c.re/N;
				h.im += c.im/N;
			} else {
				h.re += c.re;
				h.im += c.im;
			}
		}
		results[n] = h;
	}
	return results;
}


// ### Reading .wav Data ###


#pragma pack(push, 1)
typedef struct {
    char chunkID[4];       // "RIFF"
    uint32_t chunkSize;    // Size of the file
    char format[4];        // "WAVE"
    char subchunk1ID[4];   // "fmt "
    uint32_t subchunk1Size;// Size of the fmt chunk
    uint16_t audioFormat;  // Audio format (1 = PCM)
    uint16_t numChannels;  // Number of channels
    uint32_t sampleRate;   // Sampling rate (e.g., 44100)
    uint32_t byteRate;     // byte rate = sampleRate * numChannels * bitsPerSample/8
    uint16_t blockAlign;   // block align = numChannels * bitsPerSample/8
    uint16_t bitsPerSample;// Bits per sample (8, 16, etc.)
    char subchunk2ID[4];   // "data"
    uint32_t subchunk2Size;// Number of bytes in the data
} WAVHeader;
#pragma pack(pop)

void printWAVHeader(WAVHeader header) {
    printf("Chunk ID: %.4s\n", header.chunkID);
    printf("Chunk Size: %u\n", header.chunkSize);
    printf("Format: %.4s\n", header.format);
    printf("Subchunk1 ID: %.4s\n", header.subchunk1ID);
    printf("Subchunk1 Size: %u\n", header.subchunk1Size);
    printf("Audio Format: %u\n", header.audioFormat);
    printf("Number of Channels: %u\n", header.numChannels);
    printf("Sample Rate: %u\n", header.sampleRate);
    printf("Byte Rate: %u\n", header.byteRate);
    printf("Block Align: %u\n", header.blockAlign);
    printf("Bits Per Sample: %u\n", header.bitsPerSample);
    printf("Subchunk2 ID: %.4s\n", header.subchunk2ID);
    printf("Subchunk2 Size: %u\n", header.subchunk2Size);
}


// ### main ###


int main(int argc, char *argv[]) {
    
    if (argc < 2) {
        printf("Expected executable in form: %s <filename.wav>\n", argv[0]);
        return 1;
    }

    FILE *file = fopen(argv[1], "rb");
    if (!file) {
        perror("Could not open file :(");
        return 1;
    }

    WAVHeader header;

    if (fread(&header, sizeof(WAVHeader), 1, file) != 1) {
        perror("Could not read file header");
        fclose(file);
        return 1;
    }

    // Check file is RIFF or WAVE
    if (strncmp(header.chunkID, "RIFF", 4) != 0 || strncmp(header.format, "WAVE", 4) != 0) {
        printf("Not a valid .wav file.\n");
        fclose(file);
        return 1;
    }

    // Check if the subchunk1ID is "fmt "
    if (strncmp(header.subchunk1ID, "fmt ", 4) != 0) {
        printf("Error: Expected 'fmt ' subchunk not found.\n");
        fclose(file);
        return 1;
    }
    //Print .wav header for initial bugfixing
    printWAVHeader(header);

    fseek(file, sizeof(WAVHeader), SEEK_SET);
    int oneChannelLength = header.subchunk2Size * 8 / header.bitsPerSample / 2;
    int numRows = oneChannelLength * FPS / (header.sampleRate);
    int numCols = oneChannelLength / numRows / s_factor;
    size_t rows_size = numRows * sizeof(cnum *);

    cnum **leftData = (cnum**)malloc(rows_size);
    cnum **rightData = (cnum**)malloc(rows_size);

    if (leftData == NULL || rightData == NULL) {
        perror("Failed to allocate memory to audio data arrays");
        return 1;
    }

    for (int row = 0; row < numRows; row++) {
        leftData[row] = (cnum *)malloc(numCols * sizeof(cnum));
        rightData[row] = (cnum *)malloc(numCols * sizeof(cnum));
    }
    printf("\n\n%d rows, %d cols.\n\n",numRows,numCols);
    int16_t sample;

    for (int row = 0; row < numRows; row++) {
        for (int col = 0; col < numCols; col++) {
            for (int i=0; i<s_factor;i++) {
                if (fread(&sample, sizeof(int16_t), 1, file) != 1) {
                    perror("Error reading from left channel");
                    fclose(file);
                    return 1;
                }
                leftData[row][col].re += (float)sample;

                if (fread(&sample, sizeof(int16_t), 1, file) != 1) {
                    perror("Error reading from right channel");
                    fclose(file);
                    return 1;
                }
                rightData[row][col].re += (float)sample;
            }
            leftData[row][col].re /= s_factor;
            rightData[row][col].re /= s_factor;
            leftData[row][col].im = 0.0f;
            rightData[row][col].im = 0.0f;
        }
    }

    fclose(file);

    /*Output Test for now smoothened data
    float rowTime;
    float colTime;

    for (int row = 0; row < 1; row++) {
        rowTime = row/FPS;
        for (int col = 0; col < 50; col++) {
            colTime = col/(float)header.sampleRate * s_factor;
            printf("%.5f:\t%.3f + %.3fi\t%.3f + %.3fi\n", rowTime+colTime, leftData[row][col].re, leftData[row][col].im,rightData[row][col].re, rightData[row][col].im);
        }
        printf("\n");
    }
    */
    
    // ### DFT time :p ###


    FILE *file2 = fopen("wave_output.txt","w");
    if (file2 == NULL) {
    perror("Failed to open file");
    return 1;
    } 

    fprintf(file2,"time");

    for (int col=0;col<numCols/2;col++) {
        fprintf(file2,",%d",col*44100/numCols/s_factor);
    }

    fprintf(file2,"\n");


    for (int row = 0; row < numRows; row++) {
        cnum *leftDFT = transform(DFT,leftData[row],numCols);
        cnum *rightDFT = transform(DFT,rightData[row],numCols);
        fprintf(file2,"%.5f",(float)row/FPS);

        for (int i = 0; i<numCols/2; i++) {
            fprintf(file2,",%.3f",(mag(leftDFT[i]) + mag(rightDFT[i])) / 2);
        }
        fprintf(file2,"\n");
        /*Output Test
        for (int i=0;i<numCols;i++) {
        printf("%d: %.3f\n",i*44100/numCols/s_factor,mags[row][i]);
        }
        */

        free(leftDFT);
        free(rightDFT);
    }

    fclose(file2);

    for (int row = 0; row < numRows; row++) {
        free(leftData[row]);
        free(rightData[row]);
    }
    
    
    free(leftData);
    free(rightData);

    printf("Done!\n");
    return 0;
}
