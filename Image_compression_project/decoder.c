#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "myBMPLib.h"
#include "my2DDCTLib.h"



void usage(FILE *fp)
{
    fprintf(fp, "decoder output.bin output.bmp\n");
    return;
}

int demo0(char *bmp_fn, char *fn_R, char *fn_G, char *fn_B, char *fn_dim) {
    BMP_file bmp;

    FILE *file = fopen(fn_dim, "r");
    int file_size;
    if (file == NULL) {
        fprintf(stderr, "無法打開文件 %s\n", fn_dim);
        return 1;
    }

    char buffer[256];
    if (fgets(buffer, sizeof(buffer), file) != NULL) {
        // 將字符串轉換為整數
        file_size = atoi(buffer);
        fclose(file);
    }
    FILE* outfile = fopen(bmp_fn,"wb");

    text_save_BMP(&bmp,outfile, fn_R, fn_G, fn_B, file_size);

    fclose(outfile);
    return 1;
}



int demo1(char *Qbmp,char *Qt_Y,char *Qt_Cb,char *Qt_Cr,char *fn_dim,char *qF_Y,char *qF_Cb,char *qF_Cr,char *eF_Y,char *eF_Cb,char *eF_Cr){
    BMP_file bmp;
    YCbCr value;
    my2DDCT dct;
    FILE *file = fopen(fn_dim, "r");
    ckeck_size(&bmp,&dct, &value,fn_dim);//return file_size
    cut8_8_initial(&bmp,&value);
    my2DDCT_initial(&dct, &value);
    input(&value,Qt_Y,Qt_Cb,Qt_Cr);
    input_quantize(&value, &dct,qF_Y,qF_Cb,qF_Cr,eF_Y,eF_Cb,eF_Cr);
    I_quantize(&value,&dct);
    IDCT(&value,&dct);
    initial_RGB_YCbCr(&bmp,&value);
    I_split(&bmp,&value);
    YCbCr_to_RGB(&bmp);
    RGB_to_file(&bmp,Qbmp ,ckeck_size(&bmp,&dct, &value,fn_dim));
    my2DDCT_free(&dct,&value);
    cut8_8_free(&value);
    return 1;
}

int main(int argc, char **argv)
{
    int option = 0;// 0 for demo0, 1 for demo1, etc
    BMP_file bmp; // BMP file structure

    option = atoi(argv[1]);

    switch(option)
    {
        case 0:
        demo0(argv[2], argv[3], argv[4], argv[5], argv[6]);
        break;

        case 1:
        demo1(argv[2], argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9],argv[10],argv[11],argv[12]);
        break;       
        default:
        break;
        // default statements
    }
    return 0; // 0: exit normally 1: exit with error
}