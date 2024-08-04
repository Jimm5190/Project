#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#define PI 3.14159265359

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

void writeWavHeader(FILE* file, int fs, int m, double T);
double wavefunction(char *wavetype, size_t n, double A, int f, double Ts);

int main(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s fs m wavetype f A T 1> fn.wav 2> sqnr.txt\n", argv[0]);
        return 1;
    }

    int fs = atoi(argv[1]);
    int m = atoi(argv[2]);
    char *wavetype = argv[3];
    int f = atoi(argv[4]);
    double A = atof(argv[5]);
    double T = atof(argv[6]);

#ifdef _WIN32
    // Windows-specific code
    if (_setmode(fileno(stdout), O_BINARY) == -1) {
        perror("Error setting stdout to binary mode");
        return 1;
    }
#endif

    writeWavHeader(stdout, fs, m, T) ;

    double Ts = 1.0/fs ; //因預設時間長度為T，改由Ts代表取樣週期
    size_t n ; 
    size_t N = (size_t)( 1 * T * fs );// length of sin wave (sample)
    unsigned char* x8 = NULL;
    short* x16 = NULL ;
    int32_t* x32 = NULL ;

    double maxAmplitude = 0;
    double tmp ; 
    double sqnr =0.0;
    double sum_x = 0 , sum_e = 0 ;
    double power_x = 0 , power_e = 0 ;
    double* e = (double*)malloc(sizeof(double) * N); 
   
    //建立 x 數组並進行動態分配
    if(m == 8){
        x8 = (unsigned char*)malloc(sizeof(unsigned char) * N);
    }else if(m == 16){
        x16 = (short*)malloc(sizeof(short) * N);
    }else if(m == 32){
        x32 = (int32_t*)malloc(sizeof(int32_t) * N);
    }    
    //取樣      
    maxAmplitude = pow(2.0 , m-1)-1;
    for( n = 0; n < N; n ++ ) {
        //產生該wave第n個取樣點的amplitude
        tmp = maxAmplitude * wavefunction(wavetype,n,A,f,Ts) ;
        
        //quantization  & calculate SQNR
        if(m == 8){
            tmp +=128.0 ;//sample size = 8bits時範圍0 ~ 255
            x8[n] = (unsigned char)floor(tmp + 0.5);
            e[n] = tmp - (double)x8[n];
        } else if(m == 16){
            x16[n] = (short)floor(tmp + 0.5);
            e[n] = tmp - (double)x16[n] ;
        } else if(m == 32){
            x32[n] = (int32_t)floor(tmp + 0.5);
            e[n] = tmp - (double)x32[n] ;
        }
        sum_x += tmp * tmp ;
        sum_e += e[n] * e[n] ;
    }

    power_x = sum_x/N ;
    power_e = sum_e/N ;
    sqnr = 10 * log10(power_x/power_e) ;
    sqnr = round(sqnr * pow(10 , 15)) / pow(10 , 15) ;

    fprintf(stderr, "SQNR : %.15f", sqnr);

    fflush(stderr);

    if(m == 8){
        fwrite(x8, sizeof(unsigned char), N, stdout);
        free(x8) ;
    } else if(m == 16){
        fwrite(x16, sizeof(short), N, stdout);    
        free(x16) ;
    }else if(m == 32){
        fwrite(x32, sizeof(int32_t), N, stdout);
        free(x32);
    }

    free(e) ;
    fclose(stdout);
    fclose(stderr);

    return 0;
}

void writeWavHeader(FILE* file, int fs, int m, double T) {
    uint8_t header[44] = {
        'R', 'I', 'F', 'F',
        0, 0, 0, 0,
        'W', 'A', 'V', 'E',
        'f', 'm', 't', ' ',
        16, 0, 0, 0,
        1, 0,
        1, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0,
        0, 0,
        'd', 'a', 't', 'a',
        0, 0, 0, 0
    };

    int byteRate = (fs * 1 * m) / 8;
    int subChunk2Size = byteRate * T;
    int ChunkSize = subChunk2Size + 36;

    header[4] = ChunkSize & 0xFF;
    header[5] = (ChunkSize >> 8) & 0xFF;
    header[6] = (ChunkSize >> 16) & 0xFF;
    header[7] = (ChunkSize >> 24) & 0xFF;
    header[24] = fs & 0xFF;
    header[25] = (fs >> 8) & 0xFF;
    header[26] = (fs >> 16) & 0xFF;
    header[27] = (fs >> 24) & 0xFF;
    header[28] = byteRate & 0xFF;
    header[29] = (byteRate >> 8) & 0xFF;
    header[30] = (byteRate >> 16) & 0xFF;
    header[31] = (byteRate >> 24) & 0xFF;
    header[32] = (1 * m) / 8;
    header[34] = m;
    header[40] = subChunk2Size & 0xFF;
    header[41] = (subChunk2Size >> 8) & 0xFF;
    header[42] = (subChunk2Size >> 16) & 0xFF;
    header[43] = (subChunk2Size >> 24) & 0xFF;

    fwrite(header, sizeof(uint8_t), 44, file);
}

double wavefunction(char *wavetype, size_t n, double A, int f, double Ts){
    double tmp ;
    if(strcmp(wavetype, "sine") == 0){ 
        tmp = A * sin(2 * PI * f * n * Ts);
    } else if(strcmp(wavetype, "square") == 0){ 
        if(fmod(n * Ts * f, 1.0) < 0.5){
            tmp = A;
        } else{
            tmp = -1 * A;
        }       
    } else if(strcmp(wavetype, "triangle") == 0){
        if(fmod(n * Ts * f, 1.0) <= 0.25 ) {
            tmp = A * (4 * fmod(n * Ts * f, 1.0));
        } else if(fmod(n * Ts * f, 1.0) <= 0.75) {
            tmp = A * (-4 * (fmod(n * Ts * f, 1.0) - 0.25) + 1);
        } else {
            tmp = A * (4 * (fmod(n * Ts * f, 1.0) - 1.0));
        }
    } else if(strcmp(wavetype, "sawtooth") == 0){
        if(fmod(n * Ts * f, 1.0) < 1 ) {
            tmp = A * (2 * f * n * Ts - 2 * (int)(f * n * Ts) - 1);
        } else{
            tmp = A;
        }
    }
    return tmp;
}
