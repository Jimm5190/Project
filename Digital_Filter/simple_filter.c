#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>

#define FS 48000.0f
#define FL 1500.0f
#define FH 3500.0f
#define PI 3.141592653589793f
#define N_total 1200


typedef struct _wav {
	int fs;
	char header[44];
	size_t length;
	short *LChannel;
	short *RChannel;
} wav;

int wav_read_fn(char *fn, wav *p_wav)
{
	//char header[44];
	short temp = 0;
	size_t i = 0;

	FILE *fp = fopen(fn, "rb");
	if(fp==NULL) {
		fprintf(stderr, "cannot read %s\n", fn);
		return 0;
	}
	fread(p_wav->header, sizeof(char), 44, fp);
	while( !feof(fp) ) {
		fread(&temp, sizeof(short), 1, fp);
		i++;
	}
	p_wav->length = i / 2;
	p_wav->LChannel = (short *) calloc(p_wav->length, sizeof(short));
	if( p_wav->LChannel==NULL ) {
		fprintf(stderr, "cannot allocate memory for LChannel in wav_read_fn\n");
		fclose(fp);
		return 0;
	}
	p_wav->RChannel = (short *) calloc(p_wav->length, sizeof(short));
	if( p_wav->RChannel==NULL ) {
		fprintf(stderr, "cannot allocate memory for RChannel in wav_read_fn\n");
		fclose(fp);
		return 0;
	}
	fseek(fp, 44, SEEK_SET);
	for(i=0;i<p_wav->length;i++) {
		fread(p_wav->LChannel+i, sizeof(short), 1, fp);
		fread(p_wav->RChannel+i, sizeof(short), 1, fp);
	}
	fclose(fp);
	return 1;
}

int wav_save_fn(char *fn, wav *p_wav)
{
	FILE *fp = fopen(fn, "wb");
	size_t i;
	if(fp==NULL) {
		fprintf(stderr, "cannot save %s\n", fn);
		return 0;
	}
	fwrite(p_wav->header, sizeof(char), 44, fp);
	for(i=0;i<p_wav->length;i++) {
		fwrite(p_wav->LChannel+i, sizeof(short), 1, fp);
		fwrite(p_wav->RChannel+i, sizeof(short), 1, fp);
	}
	fclose(fp);
	return 1;
}

int wav_init(size_t length, wav *p_wav)
{
	p_wav->length = length;
	p_wav->LChannel = (short *) calloc(p_wav->length, sizeof(short));
	if( p_wav->LChannel==NULL ) {
		fprintf(stderr, "cannot allocate memory for LChannel in wav_read_fn\n");
		return 0;
	}
	p_wav->RChannel = (short *) calloc(p_wav->length, sizeof(short));
	if( p_wav->RChannel==NULL ) {
		fprintf(stderr, "cannot allocate memory for RChannel in wav_read_fn\n");
		return 0;
	}
	return 1;
}

void wav_free(wav *p_wav)
{
	free(p_wav->LChannel);
	free(p_wav->RChannel);
}

/* hamming: for n=0,1,2,...N, length of N+1 */
float hamming(int N, int n)
{
	return 0.54 - 0.46 * cosf(2*PI*((float)(n))/((float)N));
}

/* low-pass filter coef: n=0,1,2...,2M */
float low_pass(int m, int n)
{
	float wc = 2*PI*FL/FS;
	if(n==m) {// L'Hopital's Rule
		return wc/PI;
	}
	else {
		return sinf(wc*((float)(n-m)))/PI/((float)(n-m)) * hamming(2*m+1, n);
	}
}

float high_pass(int m, int n)
{
	float wc = 2*PI*FL/FS;
	if(n==m) {// L'Hopital's Rule
		return 2.0*(1.0-(wc/PI));
	}
	else {
		return -2.0*sinf(wc*((float)(n-m)))/PI/((float)(n-m)) * hamming(2*m+1, n);
	}
}

float band_pass(int m, int n)
{
	float wh = 2*PI*FH/FS;
    float wl = 2*PI*FL/FS;
	if(n==m) {// L'Hopital's Rule
		return 2.0*(wh/PI - wl/PI);
	}
	else {
		return 2.0*(sinf(wh*((float)(n-m)))-sinf(wl*((float)(n-m))))/PI/((float)(n-m)) * hamming(2*m+1, n);
	}
}


float band_stop_pass(int m, int n)
{
	float wh = 2*PI*FH/FS;
    float wl = 2*PI*FL/FS;
	if(n==m) {// L'Hopital's Rule
		return 1.0*(1.0- (wh/PI - wl/PI));
	}
	else {
		return 1.0*(-(sinf(wh*((float)(n-m)))-sinf(wl*((float)(n-m))))/PI/((float)(n-m)) )* hamming(2*m+1, n);
	}
}



int main(int argc, char *argv[])
{
	wav wavin;
	wav wavout;
	int M = atoi(argv[1]);
	char hL_txt[1024] ;
	strncpy(hL_txt, argv[2], sizeof(hL_txt) - 1);
    hL_txt[sizeof(hL_txt) - 1] = '\0';
	char hR_txt[1024] ;
	strncpy(hR_txt, argv[3], sizeof(hR_txt) - 1);
    hR_txt[sizeof(hR_txt) - 1] = '\0';
	char YL_txt[1024] ;
	strncpy(YL_txt, argv[4], sizeof(YL_txt) - 1);
    YL_txt[sizeof(YL_txt) - 1] = '\0';
	char YR_txt[1024] ;
	strncpy(YR_txt, argv[5], sizeof(YR_txt) - 1);
    YR_txt[sizeof(YR_txt) - 1] = '\0';
	char fn_in[1024] ;
	strncpy(fn_in, argv[6], sizeof(fn_in) - 1);
    fn_in[sizeof(fn_in) - 1] = '\0';
	char fn_out[1024] ;
	strncpy(fn_out, argv[7], sizeof(fn_out) - 1);
    fn_out[sizeof(fn_out) - 1] = '\0';


	float h_L[2050] = {0};
    float h_R[2050] = {0};
	double exp_cos[N_total] ={0};
	double exp_sin[N_total] ={0};
	double y_L[N_total] = {0};
	double y_R[N_total] = {0};
	double YL_re[N_total/2+1] ={0};
	double YL_im[N_total/2+1] ={0};
	double YR_re[N_total/2+1] ={0};
	double YR_im[N_total/2+1] ={0};
	int n = 0;
	float y = 0;
	double tmp_L = 0,tmp_R = 0;
	double magnitude = 0;
	int k;


	// read wav
	if( wav_read_fn(fn_in, &wavin) == 0 ) {
		fprintf(stderr, "cannot read wav file %s\n", fn_in);
		exit(1);
	}


	// construct low-pass filter
	for(n=0;n<(2*M+1);n++) {
		h_L[n] = band_pass(M, n);
        h_R[n] = band_stop_pass(M, n);
	}

    FILE* hL = fopen(hL_txt , "w") ;
	FILE* hR = fopen(hR_txt , "w") ;

	for(n=0;n<(2*M+1);n++) {
		fprintf(hL, "%e\n", h_L[n]);
		fprintf(hR, "%e\n", h_R[n]);
	}

	fclose(hL);
	fclose(hR);
    

	// filtering (convolution)
	if( wav_init(wavin.length, &wavout)==0 ) {
		exit(1);
	}

	for(n=0;n<wavin.length;n++) {
		y = 0;
		for(k=0;k<(2*M+1);k++) {
			if( (n-k)>=0 )
				y = y + h_L[k] * ((float)(wavin.LChannel[n-k]));
		}
		if(n>=962880 && n<964080){//20.060 秒到 20.085 秒之間的 1200 個 sample 點
			y_L[n-962880] = (float)(roundf(y)) ;
		}
		wavout.LChannel[n] = (short)(roundf(y));

		y = 0;
		for(k=0;k<(2*M+1);k++) {
			if( (n-k)>=0 )
				y = y + h_R[k] * ((float)(wavin.RChannel[n-k]));
		}
		if(n>=962880 && n<964080){
			y_R[n-962880] = (float)(roundf(y))  ;
		}
		wavout.RChannel[n] = (short)(roundf(y));
	}
	memcpy(wavout.header, wavin.header, 44);
	//DFT
	for(n=0;n<N_total;n++){
	    exp_cos[n] = cosf(2.0 * PI / N_total  * n); //先計算exp後的值減少計算複雜度
        exp_sin[n] = sinf(2.0 * PI / N_total  * n);
    }
	for(k=0;k<(N_total)/2+1;k++) {
	    YL_re[k] = 0.0;
	    YL_im[k] = 0.0;
	    for(n=0;n<N_total;n++) {	    
			tmp_L = y_L[n] * exp_cos[(n*k) % N_total] ;
			YL_re[k] += tmp_L;
			tmp_L =  y_L[n] * -exp_sin[(n*k) % N_total];
			YL_im[k] += tmp_L;
	    	}
        YL_re[k] /= N_total; //normalize
		YL_im[k] /= N_total; //normalize
    }

	for(k=0;k<(N_total)/2+1;k++) {
	    YR_re[k] = 0;
	    YR_im[k] = 0;
	    for(n=0;n<N_total;n++) {	    
		tmp_R = y_R[n] * exp_cos[(n*k)% N_total] ;
		YR_re[k] += tmp_R;
		tmp_R = y_R[n] * -exp_sin[(n*k)% N_total];
		YR_im[k] += tmp_R;
	    }
        YR_re[k] /= N_total; //normalize
		YR_im[k] /= N_total; //normalize
    }

	FILE* YL = fopen(YL_txt , "w") ;
	FILE* YR = fopen(YR_txt , "w") ;
    //轉換成dB單位
    for(k = 0; k < (N_total)/2+1; k ++ ) {  
	    magnitude = 20.0 * log10( YL_re[k] * YL_re[k] + YL_im[k] * YL_im[k]) ;
	    if(magnitude < 0)
			magnitude = 0.0 ;
		
	    fprintf(YL, "%e\n", magnitude);
	    fflush(YL);
	}
	for(k = 0; k <(N_total)/2+1; k ++ ) {  
	    magnitude = 20.0 * log10(YR_re[k] * YR_re[k] + YR_im[k] * YR_im[k]);
	    if(magnitude < 0)
			magnitude = 0.0 ;
		
	    fprintf(YR, "%e\n", magnitude);
	    fflush(YR);
	}

	fclose(YL);
	fclose(YR);

	// save wav
	if( wav_save_fn(fn_out, &wavout)==0) {
		fprintf(stderr, "cannot save %s\n", fn_out);
		exit(1);

	}
	wav_free(&wavin);
	wav_free(&wavout);
}



