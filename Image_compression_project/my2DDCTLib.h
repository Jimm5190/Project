#ifndef INCLUDE_my2DDCTLib
#define INCLUDE_my2DDCTLib
#include "myBMPLib.h"

typedef struct _my2DDCT {
	size_t H; // height
	size_t W; // width
	float ****basis;
	float ****Y_basis;
	float ****Cb_basis;
	float ****Cr_basis;
	int ****Y;//8*8 F[u,v]quantizate
	int ****Cb;//8*8 F[u,v]quantizate
	int ****Cr;//8*8 F[u,v]quantizate
	float ****eY;
	float ****eCb;
	float ****eCr;
	//float basis[8][8][8][8];
} my2DDCT;

typedef struct _YCbCr{
	size_t count_H; // height
	size_t count_W; // width
	float ****Y;//8*8 f[r,c]
	float ****Cb;//8*8 f[r,c]
	float ****Cr;//8*8 f[r,c]
	int **Y_qt;
	int **Cb_qt;
	int **Cr_qt;
	int ***Y_zigzag;
	int ***Cb_zigzag;
	int ***Cr_zigzag;
	int ***Y_rle;
	int ***Cb_rle;
	int ***Cr_rle;
	int *Y_rle_code;//merge
	int *Cb_rle_code;//merge
	int *Cr_rle_code;//merge
	int Y_block_size ;
	int Cb_block_size ;
	int Cr_block_size ;
}YCbCr;

typedef struct _RLE
{
	int *Y_num ;
	int *Y_count ;
	int *Cb_num ;
	int *Cb_count ;
	int *Cr_num ;
	int *Cr_count ;
	int Y_len;
	int Cb_len;
	int Cr_len;
}RLE;
int RGB_to_YCbCr(BMP_file *p_bmp);
/* decode時會沒有bmp可以assign給count_H，特地獨立出來 */
int encode_size(BMP_file *p_bmp,YCbCr *p_value);
/* 分割成每塊8*8 */
int split(BMP_file *p_bmp,YCbCr *p_value);
/* DCT */
int DCT(my2DDCT *p_dct,YCbCr *p_value);
/* 以 ascii 儲存的channel  的 quantization table & 將值存到cut8_8_qt struct中*/
int quantization_table(YCbCr *p_value);
/* calculate SQNR*/
int SQNR(my2DDCT *p_dct,YCbCr *p_value);
/* output_quantize_table*/
int output_quantize_table(char *Qt_Y,char *Qt_Cb,char *Qt_Cr,YCbCr *p_value);
/* DCT後的值quantization */
int quantize(my2DDCT *p_dct, YCbCr *p_value);
/* DPCM DC*/
int DPCM(my2DDCT *p_dct,YCbCr *p_value);
/* zigzag */
int zigzag(my2DDCT *p_dct,YCbCr *p_value);
/* Run Length Encoding (RLE) */
int RLE1(YCbCr *p_value);
/* RLE2 */
int RLE2(YCbCr *p_value, RLE *p_rle);

//initial
/* initialize 2D-DCT bases */
int my2DDCT_initial(my2DDCT *p_dct,YCbCr *p_value);
/* initialize cut8_8  */
int cut8_8_initial(BMP_file *p_bmp,YCbCr *p_value);

//free
/* free 2D-DCT bases */
int my2DDCT_free(my2DDCT *p_dct,YCbCr *p_value);
/* free cut8_8_qt */
int cut8_8_free(YCbCr *p_value);

//output
/* output_quantize */
int output_quantize(YCbCr *p_value,my2DDCT *p_dct,char *qF_Y,char* qF_Cb,char* qF_Cr);
/* output_quantize_error */
int output_quantize_error(YCbCr *p_value, my2DDCT *p_dct,char *eF_Y,char* eF_Cb,char* eF_Cr);
/* output rle*/
int output_rle(BMP_file *p_bmp,YCbCr *p_value,char *rle_code);

//decode
/* input file*/
int input(YCbCr *p_value,char *Qt_Y,char *Qt_Cb,char *Qt_Cr);
int input_quantize(YCbCr *p_value, my2DDCT *p_dct,char *qF_Y,char* qF_Cb,char* qF_Cr,char *eF_Y,char *eF_Cb,char *eF_Cr);
int I_quantize(YCbCr *p_value,my2DDCT *p_dct);
int IDCT(YCbCr *p_value,my2DDCT *p_dct);
int I_split(BMP_file *p_bmp,YCbCr *p_value);
int YCbCr_to_RGB(BMP_file *p_bmp);
int RGB_to_file(BMP_file *p_bmp,char *bmp_fn ,int file_size);
int initial_RGB_YCbCr(BMP_file *p_bmp,YCbCr *p_value);
int ckeck_size(BMP_file *p_bmp,my2DDCT *p_dct, YCbCr *p_value,char*fn_dim);
#endif