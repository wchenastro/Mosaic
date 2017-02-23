#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>
#include <math_constants.h>
#include <vector_types.h>
//#include <cuComplex.h>

//#define speedOfLight 299792458

__global__ void waveNumberMatrix(double* waveNumbers, double* horizontal,
        double* frequencies, int beamNumber, int freqNumber) {

    const long speedOfLight = 299792458;
    int beamIdx = threadIdx.x;
    int waveLengthIdx = blockIdx.x;
    if(beamIdx > beamNumber or waveLengthIdx > freqNumber) return;

    double altitude = horizontal[beamIdx*2];
    double azimuth = horizontal[beamIdx*2 + 1];
    double theta = CUDART_PI/2. - altitude;
    double phi = CUDART_PI/2.- azimuth;
    long freqChannel = waveLengthIdx*beamNumber*3;

    double wavenumber = (-1) * 2 * CUDART_PI * frequencies[waveLengthIdx] / speedOfLight;

    waveNumbers[freqChannel +       0        + beamIdx] = sin(theta) * cos(phi) * wavenumber;
    waveNumbers[freqChannel + beamNumber * 1 + beamIdx] = sin(theta) * sin(phi) * wavenumber;
    waveNumbers[freqChannel + beamNumber * 2 + beamIdx] = cos(theta) * wavenumber;

    return;
}


__global__ void convertTo8Bits(int8_t* convertedWeights, double* weights, int freqChannels, float scale){

    volatile int matrixRightWidth = blockDim.x;
    volatile int matrixLeftRowIdx = blockIdx.x;
    volatile int matrixRightColIdx = threadIdx.x;
    volatile unsigned int weightsIndex = blockIdx.y * gridDim.x * matrixRightWidth *2 +
        matrixLeftRowIdx * matrixRightWidth * 2 + matrixRightColIdx * 2;
    volatile double weightReal = weights[weightsIndex];
    volatile double weightImag = weights[weightsIndex + 1];
    volatile double absolteWeight = sqrt(weightReal*weightReal + weightImag*weightImag);
    convertedWeights[weightsIndex] = (int8_t) scale * (weightReal / absolteWeight);
    convertedWeights[weightsIndex + 1] = (int8_t) scale * (weightImag / absolteWeight);

    return;
}

__global__ void weightMatrix(double* weights, double* matrixLeft, double* matrixRight,
        int matrixLeftWidth, int matrixRightWidth, int freqChannels) {

    double vectorSum = 0.0;
    int cursor = 0;
    double complexNumber[2];

    volatile int matrixLeftRowIdx = blockIdx.x;
    volatile int matrixRightColIdx = threadIdx.x;
    volatile int subMatrixRightStarter = blockIdx.y * matrixLeftWidth * matrixRightWidth;
    volatile int matrixLeftRowStarter = matrixLeftRowIdx * matrixLeftWidth;

    for(cursor = 0; cursor < matrixLeftWidth; cursor++){
        vectorSum += matrixLeft[matrixLeftRowStarter + cursor] *
            matrixRight[subMatrixRightStarter +
            matrixRightWidth * cursor + matrixRightColIdx];
    }
    sincos(-vectorSum, complexNumber, complexNumber+1);
    volatile unsigned int weightsIndex = blockIdx.y * gridDim.x * matrixRightWidth *2 +
        matrixLeftRowIdx * matrixRightWidth * 2 + matrixRightColIdx * 2;
    weights[weightsIndex] = complexNumber[1];
    weights[weightsIndex + 1] = complexNumber[0];

    return;
}


__global__ void weightMatrix8BitsScaled(char2* weights, double* matrixLeft, double* matrixRight,
        int matrixLeftWidth, int matrixRightWidth, int freqChannels, float scale) {

    double vectorSum = 0.0;
    int cursor = 0;
    double2 complexNumber;

    volatile int matrixLeftRowIdx = blockIdx.x;
    volatile int matrixRightColIdx = threadIdx.x;
    volatile int subMatrixRightStarter = blockIdx.y * matrixLeftWidth * matrixRightWidth;
    volatile int matrixLeftRowStarter = matrixLeftRowIdx * matrixLeftWidth;

    for(cursor = 0; cursor < matrixLeftWidth; cursor++){
        vectorSum += matrixLeft[matrixLeftRowStarter + cursor] *
            matrixRight[subMatrixRightStarter +
            matrixRightWidth * cursor + matrixRightColIdx];
    }
    sincos(-vectorSum, &complexNumber.y, &complexNumber.x);
    volatile double absolteWeight = sqrt(pow(complexNumber.x, 2) + pow(complexNumber.y, 2));
    volatile unsigned int weightsIndex = blockIdx.y * gridDim.x * matrixRightWidth +
        matrixLeftRowIdx * matrixRightWidth + matrixRightColIdx;
    weights[weightsIndex].x = (int8_t) scale * (complexNumber.x / absolteWeight);;
    weights[weightsIndex].y = (int8_t) scale * (complexNumber.y / absolteWeight);;

    return;
}

int createFrequencies(double* frequencies, double centralFreq, int channels, double bandwidth) {
    int i;
    double startFreq = centralFreq - bandwidth/2.0;
    double step = bandwidth/2.0/(channels/2);
    for(i=0;i<channels;i++) {
        frequencies[i] = startFreq + step*i;
    }
    return 0;
}

int readAntennaPosition(const char * fileName, double * positions) {

    FILE* stream = fopen(fileName, "r");

    if(stream == NULL) {
        printf("%s does not exist!\n", fileName);
        return 1;
    }

    unsigned int lineLength = 100;
    char line[lineLength];
    unsigned int index = 0;
    //double a, b, c;
    while (fgets(line, lineLength, stream)) {
        positions[index++] = atof(strtok(line, " "));
        positions[index++] = atof(strtok(NULL, " "));
        positions[index++] = atof(strtok(NULL, " "));
        //printf("%.17lf %.17lf %.17lf\n", a,b,c);
    }
    //printf("===================\n");
    return 0;
}

int readBeamPosition(const char * fileName, double * positions) {

    FILE* stream = fopen(fileName, "r");

    if(stream == NULL) {
        printf("%s does not exist!\n", fileName);
        return 1;
    }

    unsigned int lineLength = 60;
    char line[lineLength];
    int index = 0;
    //double a, b;
    while (fgets(line, lineLength, stream)) {
        positions[index++] = atof(strtok(line, " "));
        positions[index++] = atof(strtok(NULL, " "));
        //printf("%.17lf %.17lf\n", a,b);
    }
    //printf("===================\n");

    return 0;
}

int writeWeigts(const char * fileName, double * weights,
        int antennaNumber, int beamNumber, int freqChannels) {

    int ant, beam, channel;
    int beamStarter = 0;
    int channelStarter = 0;

    FILE* fp = fopen(fileName, "w");

    for(channel=0;channel<freqChannels;channel++) {
        channelStarter  = channel*antennaNumber*beamNumber*2;
        for(ant=0;ant<antennaNumber;ant++) {
            beamStarter  = ant*beamNumber*2;
            for(beam=0;beam<beamNumber;beam++) {
                fprintf(fp, "%.17lf %.17lf ",
                        weights[channelStarter + beamStarter + beam*2],
                        weights[channelStarter + beamStarter + beam*2 + 1]);
            }
            fprintf(fp, "\n");
        }
    }

    fclose(fp);

    return 0;
}

int writeWeigts8Bits(const char * fileName, char2 * weights,
        int antennaNumber, int beamNumber, int freqChannels) {

    int ant, beam, channel;
    int beamStarter = 0;
    int channelStarter = 0;

    FILE* fp = fopen(fileName, "w");

    for(channel=0;channel<freqChannels;channel++) {
        channelStarter  = channel*antennaNumber*beamNumber;
        for(ant=0;ant<antennaNumber;ant++) {
            beamStarter  = ant*beamNumber;
            for(beam=0;beam<beamNumber;beam++) {
                fprintf(fp, "%d %d ",
                        weights[channelStarter + beamStarter + beam].x,
                        weights[channelStarter + beamStarter + beam].y);
            }
            fprintf(fp, "\n");
        }
    }

    fclose(fp);

    return 0;
}


int writeWavenumber(const char * fileName, double * waveNumbers,
        int antennaNumber, int beamNumber, int freqChannels) {

    int ant, beam, channel;
    int beamStarter = 0;
    int channelStarter = 0;

    FILE* fp = fopen(fileName, "w");

    for(channel=0;channel<freqChannels;channel++) {
        channelStarter  = channel*3*beamNumber;
        for(ant=0;ant<3;ant++) {
            beamStarter  = ant*beamNumber;
            for(beam=0;beam<beamNumber;beam++) {
                fprintf(fp, "%.17lf ", waveNumbers[channelStarter + beamStarter + beam]);
            }
            fprintf(fp, "\n");
        }
    }

    fclose(fp);

    return 0;
}


unsigned int countLineOfFile(const char * fileName) {
    int c;
    unsigned int count = 0;

    FILE* stream = fopen(fileName, "r");

    if(stream == NULL) {
        printf("%s does not exist!\n", fileName);
        exit(0);
    }

    while ((c=fgetc(stream)) != EOF) {
        if ( c == '\n' )
            count++;
    }

    fclose(stream);

    return count;
}



int main(void) {

    //struct timeval transferStart, transferEnd, calcStart, calcEnd;

    char beamCoorFile[] = "beamCoor";
    char antennaCoorFile[] = "antennaCoor";
    //char weightsFile[] = "weightsCUDA";
    int freqChannels = 64;
    float scale = 127.;
    double centralFreq = 299792458/0.21;
    double bandwidth = 2.0E+6;
    double frequencies[freqChannels];
    unsigned int beamNumber = countLineOfFile(beamCoorFile);
    unsigned int antennaNumber = countLineOfFile(antennaCoorFile);
    double beamPosition[beamNumber*2];
    double antennaPosition[antennaNumber*3];
    //double waveNumbers[beamNumber*3*freqChannels];
    //double weights[freqChannels*beamNumber*2*antennaNumber];
    //long weightsSize = freqChannels*beamNumber*2*antennaNumber*sizeof(double);
    long weightsSize8Bits = freqChannels*beamNumber*antennaNumber*sizeof(char2);
    //double* weights = (double* ) malloc(weightsSize);
    char2* weights8Bits = (char2* ) malloc(weightsSize8Bits);
    readBeamPosition(beamCoorFile, beamPosition);
    readAntennaPosition(antennaCoorFile, antennaPosition);

    createFrequencies(frequencies, centralFreq, freqChannels, bandwidth);

    double *deviceBeamPosition;
    double *deviceAntennaPosition;
    //double *deviceWeights;
    char2 *deviceWeights8Bits;
    double *deviceWavenumbers;
    double *deviceFrequencies;
    cudaMalloc(&deviceBeamPosition, sizeof(beamPosition));
    cudaMalloc(&deviceAntennaPosition, sizeof(antennaPosition));
    //cudaMalloc(&deviceWeights, weightsSize);
    cudaMalloc(&deviceWavenumbers, freqChannels*beamNumber*3*sizeof(double));
    cudaMalloc(&deviceFrequencies, sizeof(frequencies));
    cudaMalloc(&deviceWeights8Bits, weightsSize8Bits);

    //gettimeofday(&transferStart, 0);
    cudaMemcpy(deviceBeamPosition, beamPosition, sizeof(beamPosition), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceAntennaPosition, antennaPosition, sizeof(antennaPosition), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceFrequencies, frequencies, sizeof(frequencies), cudaMemcpyHostToDevice);
    //cudaThreadSynchronize();
    //gettimeofday(&transferEnd, 0);

    //gettimeofday(&calcStart, 0);
    waveNumberMatrix<<<freqChannels, beamNumber>>>(deviceWavenumbers,
            deviceBeamPosition, deviceFrequencies, beamNumber, freqChannels);

    dim3 dimGrid(antennaNumber, freqChannels);
    weightMatrix8BitsScaled<<<dimGrid, beamNumber>>>(deviceWeights8Bits,
            deviceAntennaPosition, deviceWavenumbers, 3, beamNumber, freqChannels, scale);

    //convertTo8Bits<<<dimGrid, beamNumber>>>(deviceWeights8Bits, deviceWeights,
    //         freqChannels, scale);

    //cudaDeviceSynchronize();
    //gettimeofday(&calcEnd, 0);

    //cudaMemcpy(weights, deviceWeights, weightsSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(weights8Bits, deviceWeights8Bits, weightsSize8Bits, cudaMemcpyDeviceToHost);
    //cudaMemcpy(waveNumbers, deviceWavenumbers, sizeof(waveNumbers), cudaMemcpyDeviceToHost);

    //writeWeigts(weightsFile, weights, antennaNumber, beamNumber, freqChannels);
    writeWeigts8Bits("weightsCUDA8Bits", weights8Bits, antennaNumber, beamNumber, freqChannels);
    //writeWavenumber("wavenumber", waveNumbers, antennaNumber, beamNumber, freqChannels);

    cudaFree(deviceWavenumbers);
    //cudaFree(deviceWeights);
    cudaFree(deviceWeights8Bits);
    cudaFree(deviceBeamPosition);
    cudaFree(deviceAntennaPosition);
    cudaFree(deviceFrequencies);

    //free(weights);
    free(weights8Bits);

    //double transferTime = (1000000.0*(transferEnd.tv_sec-transferStart.tv_sec)
    //        + transferEnd.tv_usec-transferStart.tv_usec);

    //double calcTime = (1000000.0*(calcEnd.tv_sec-calcStart.tv_sec)
    //        + calcEnd.tv_usec-calcStart.tv_usec);

    //printf("transfer time: %lf us\ncalculation time: %lf us\n", transferTime, calcTime);

    return 0;
}
