#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "myBMPLib.h"
#include "my2DDCTLib.h"



void usage(FILE *fp)
{
    fprintf(fp, "encoder 0 input.bmp R.txt G.txt B.txt dim.txt\n");
    return;
}

int demo0(char *fn_bmp, char *fn_R, char *fn_G, char *fn_B, char *fn_dim)
{
    BMP_file bmp;
    if( BMP_file_load_fn(fn_bmp, &bmp)!= 1 ) {
        exit(1);
    }
    // do something with bmp ==> datum for R, G, B, and dim
    BMP_save_text(&bmp, fn_R, fn_G, fn_B, fn_dim);

    return 1;
}


int demo1(char *fn_bmp, char *Qt_Y,char *Qt_Cb,char *Qt_Cr,char *fn_dim,char *qF_Y,char *qF_Cb,char *qF_Cr,char *eF_Y,char *eF_Cb,char *eF_Cr)
{

    BMP_file bmp;
    YCbCr value;
    my2DDCT dct;
    if( BMP_file_load_fn(fn_bmp, &bmp)!= 1 ) {
        exit(1);
    }
    
    RGB_to_YCbCr(&bmp);
    encode_size(&bmp,&value);
    cut8_8_initial(&bmp,&value);//initial 
    my2DDCT_initial(&dct, &value);
    split(&bmp,&value);
    DCT(&dct, &value);
    quantization_table(&value);
    output_quantize_table(Qt_Y,Qt_Cb,Qt_Cr,&value);
    output_size(&bmp,fn_dim);
    quantize(&dct,&value);//quantization
    output_quantize(&value,&dct,qF_Y,qF_Cb,qF_Cr);
    output_quantize_error(&value, &dct,eF_Y,eF_Cb,eF_Cr);
    SQNR(&dct, &value);
    my2DDCT_free(&dct,&value);
    cut8_8_free(&value);
    
    return 1;
}
int demo2(char *fn_bmp, char *rle_code){
    BMP_file bmp;
    YCbCr value;
    my2DDCT dct;
    RLE rle;
    if( BMP_file_load_fn(fn_bmp, &bmp)!= 1 ) {
        exit(1);
    }
    RGB_to_YCbCr(&bmp);
    encode_size(&bmp,&value);
    cut8_8_initial(&bmp,&value);//initial 
    split(&bmp,&value);//分割成每塊8*8
    my2DDCT_initial(&dct, &value);
    DCT(&dct, &value);//initial & DCT(f[r,c]-> F[u,v])
    quantization_table(&value);
    quantize(&dct,&value);//quantization
    DPCM(&dct,&value);
    zigzag(&dct,&value);
    RLE1(&value);
    output_rle(&bmp,&value,rle_code);
    RLE2(&value,&rle);
    my2DDCT_free(&dct,&value);
    cut8_8_free(&value);
}


int main(int argc, char **argv)
{
    int option = 0;// 0 for demo0, 1 for demo1, etc
    BMP_file bmp; // BMP file structure

    if ( argc <=2 ) {
        usage(stderr);
        exit(1);
    }

    option = atoi(argv[1]);


    switch(option)
    {
        case 0:
        demo0(argv[2], argv[3], argv[4], argv[5], argv[6]);
        break;

        case 1:
        demo1(argv[2], argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9],argv[10],argv[11],argv[12]);
        break;
        case 2:
        demo2(argv[2],argv[4]);
        break;
        default:
        break;
        // default statements
    }
    return 0; // 0: exit normally 1: exit with error
}