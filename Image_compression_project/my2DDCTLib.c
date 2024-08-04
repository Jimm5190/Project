#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "my2DDCTLib.h"
#include "myBMPLib.h"
#define PI 3.14159265359

int RGB_to_YCbCr(BMP_file *p_bmp)
{
    int i,j;
    if(p_bmp->H%8!=0){
        p_bmp->f_H = p_bmp->H+(8-p_bmp->H%8);//1216
    }
    if(p_bmp->W%8!=0){
        p_bmp->f_W = p_bmp->W+(8-p_bmp->W%8);//912
    }
    p_bmp->Y = (float **) calloc(p_bmp->f_H, sizeof(float *));
    for ( i=0; i<p_bmp->f_H; i++) {
        p_bmp->Y[i] = (float *) calloc(p_bmp->f_W, sizeof(float));
    }
    p_bmp->Cb = (float **) calloc(p_bmp->f_H, sizeof(float *));
    for ( i=0; i<p_bmp->f_H; i++) {
        p_bmp->Cb[i] = (float *) calloc(p_bmp->f_W, sizeof(float));
    }
    p_bmp->Cr = (float **) calloc(p_bmp->f_H, sizeof(float *));
    for ( i=0; i<p_bmp->f_H; i++) {
        p_bmp->Cr[i] = (float *) calloc(p_bmp->f_W, sizeof(float));
    }
    
    //printf("p_bmp->H: %d\n",p_bmp->H);
    for(i=0;i<p_bmp->H;i++){
        for(j=0;j<p_bmp->W;j++){
            p_bmp->Y[i][j]=0.299*(p_bmp->R[i][j])+0.587*(p_bmp->G[i][j])+0.114*(p_bmp->B[i][j]);
            p_bmp->Cb[i][j]=-0.169*(p_bmp->R[i][j])-0.331*(p_bmp->G[i][j])+0.5*(p_bmp->B[i][j])+128;
            p_bmp->Cr[i][j]=0.5*(p_bmp->R[i][j])-0.419*(p_bmp->G[i][j])-0.081*(p_bmp->B[i][j])+128;
            
        }
    }
    /*printf("encode R\n");
    for(i=0;i<8;i++){
        for(j=0;j<8;j++){
            printf("%zu\t",p_bmp->R[i][j]);
        }printf("\n");
    }*/
    
    return 1;
}

//decode時會沒有bmp可以assign給count_H，特地獨立出來
int encode_size(BMP_file *p_bmp,YCbCr *p_value){
    int n,m;
    p_value->count_H=p_bmp->f_H/8;//152
    p_value->count_W=p_bmp->f_W/8;//114
}

int YCbCr_to_RGB(BMP_file *p_bmp)
{
    int i,j;
    //p_bmp->H = p_bmp->H-6;1210
    //p_bmp->W = p_bmp->W-4;908
    p_bmp->B = (unsigned char **) calloc(p_bmp->H, sizeof(unsigned char *));
    for ( i=0; i<p_bmp->H; i++) {
        p_bmp->B[i] = (unsigned char *) calloc(p_bmp->W, sizeof(unsigned char));
    }
    p_bmp->G = (unsigned char **) calloc(p_bmp->H, sizeof(unsigned char *));
    for ( i=0; i<p_bmp->H; i++) {
        p_bmp->G[i] = (unsigned char *) calloc(p_bmp->W, sizeof(unsigned char));
    }
    p_bmp->R = (unsigned char **) calloc(p_bmp->H, sizeof(unsigned char *));
    for ( i=0; i<p_bmp->H; i++) {
        p_bmp->R[i] = (unsigned char *) calloc(p_bmp->W, sizeof(unsigned char));
    }
    //printf("decode R\n");
    for(i=0;i<p_bmp->H;i++){
        for(j=0;j<p_bmp->W;j++){
            p_bmp->R[i][j]=(unsigned char)(p_bmp->Y[i][j])+1.402*(p_bmp->Cr[i][j]-128);
            p_bmp->G[i][j]=(unsigned char)(p_bmp->Y[i][j])-0.344136*(p_bmp->Cb[i][j]-128)-0.714136*(p_bmp->Cr[i][j]-128);
            p_bmp->B[i][j]=(unsigned char)(p_bmp->Y[i][j])+1.772*(p_bmp->Cb[i][j]-128);
            //if(i<8 && j<8) printf("%zu\t",p_bmp->R[i][j]);
        }//if(i<8 && j<8) printf("\n");
    }
    /*printf("decode Y\n");
    for(i=10;i<18;i++){
        for(j=10;j<18;j++){
            printf("%f\t",p_bmp->Y[i][j]);
        }printf("\n");
    }
    printf("decode Cb\n");
    for(i=10;i<18;i++){
        for(j=10;j<18;j++){
            printf("%f\t",p_bmp->Cb[i][j]);
        }printf("\n");
    }
    printf("decode Cr\n");
    for(i=10;i<18;i++){
        for(j=10;j<18;j++){
            printf("%f\t",p_bmp->Cr[i][j]);
        }printf("\n");
    }*/
    //修正
    for(i=0;i<p_bmp->H;i++){
        for(j=0;j<p_bmp->W;j++){
            p_bmp->R[i][j]=p_bmp->R[i][j]+1;
            p_bmp->G[i][j]=p_bmp->G[i][j]+1;
            p_bmp->B[i][j]=p_bmp->B[i][j]+1;   
        }
    }
    return 1;
}

//分割成每塊8*8
int split(BMP_file *p_bmp,YCbCr *p_value){
    int i,j,x,y;
    int H=8,W=8;
    //printf("encode split\n");
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(x=0;x<H;x++){
                for(y=0;y<W;y++){
                    p_value->Y[i][j][x][y] = p_bmp->Y[x+(i*8)][y+(j*8)];
                    p_value->Cb[i][j][x][y] = p_bmp->Cb[x+(i*8)][y+(j*8)];
                    p_value->Cr[i][j][x][y] = p_bmp->Cr[x+(i*8)][y+(j*8)];
                    //if(i==100 && j==100) printf("%f\t",p_bmp->Cb[x][y]); //DEBUG
                }
                //if(i==100 && j==100 ) printf("\n");
            }
        }
    }
    return 1;
}

int I_split(BMP_file *p_bmp,YCbCr *p_value){
    int i,j,x,y;
    int H=8,W=8;
    //printf("decode split\n");
    p_value->count_H=p_bmp->f_H/8;//152
    p_value->count_W=p_bmp->f_W/8;//114
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(x=0;x<H;x++){
                for(y=0;y<W;y++){
                    p_bmp->Y[x+(i*8)][y+(j*8)]=p_value->Y[i][j][x][y];
                    p_bmp->Cb[x+(i*8)][y+(j*8)]=p_value->Cb[i][j][x][y];
                    p_bmp->Cr[x+(i*8)][y+(j*8)]=p_value->Cr[i][j][x][y];
                    //if(i==1 && j==1 ) printf("%f\t",p_value->Y[i][j][x][y]); //DEBUG
                }
                //if(i==1 && j==1) printf("\n");
            }
        }
    }
    /*printf("decode Y\n");
    for(i=0;i<9;i++){
        for(j=0;j<9;j++){
            printf("%f\t",p_bmp->Y[i][j]);
        }printf("\n");
    }*/
    return 1;
}

int initial_RGB_YCbCr(BMP_file *p_bmp,YCbCr *p_value){
    int i,j;

    //printf("\np_bmp->f_H: %d\n",p_bmp->f_H);
    p_bmp->Y = (float **) calloc(p_bmp->f_H, sizeof(float *));
    for ( i=0; i<p_bmp->f_H; i++) {
        p_bmp->Y[i] = (float *) calloc(p_bmp->f_W, sizeof(float));
    }
    p_bmp->Cb = (float **) calloc(p_bmp->f_H, sizeof(float *));
    for ( i=0; i<p_bmp->f_H; i++) {
        p_bmp->Cb[i] = (float *) calloc(p_bmp->f_W, sizeof(float));
    }
    p_bmp->Cr = (float **) calloc(p_bmp->f_H, sizeof(float *));
    for ( i=0; i<p_bmp->f_H; i++) {
        p_bmp->Cr[i] = (float *) calloc(p_bmp->f_W, sizeof(float));
    }
    return 1;
}
//DCT beta
double beta(int k){
    if(k==0)
        return 1.0/sqrt(2);
    else
        return 1.0;
}
//DCT
int DCT(my2DDCT *p_dct,YCbCr *p_value){
    //f[r,c] -> F[u,v]
    int i,j,u,v,r,c;
    int H=8,W=8;

    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<p_dct->H;u++){
                for(v=0;v<p_dct->W;v++){//p_dct->W=8;
                    for(r=0;r<p_dct->H;r++){
                        for(c=0;c<p_dct->W;c++){
                            // 計算 DCT 基函数
                            p_dct->basis[u][v][r][c]=beta(u)*beta(v)*cos(PI*u*(2.0*r+1.0)/2.0/p_dct->H)*cos(PI*v*(2.0*c+1.0)/2.0/p_dct->W);
                            // 計算 Y、Cb、Cr 基函数的 DCT 變換
                            p_dct->Y_basis[i][j][u][v]+=(p_value->Y[i][j][r][c])*p_dct->basis[u][v][r][c];
                            p_dct->Cb_basis[i][j][u][v]+=(p_value->Cb[i][j][r][c])*p_dct->basis[u][v][r][c];
                            p_dct->Cr_basis[i][j][u][v]+=(p_value->Cr[i][j][r][c])*p_dct->basis[u][v][r][c];
                
                        }
                    }
                    p_dct->Y_basis[i][j][u][v] = p_dct->Y_basis[i][j][u][v]*sqrt(4.0/H/W);
                    p_dct->Cb_basis[i][j][u][v] = p_dct->Cb_basis[i][j][u][v]*sqrt(4.0/H/W);
                    p_dct->Cr_basis[i][j][u][v] = p_dct->Cr_basis[i][j][u][v]*sqrt(4.0/H/W);
                }
            }
        }
    }
    /*printf("encode DCT\n");
    for(r=0;r<8;r++){
        for(c=0;c<8;c++){
            printf("%f\t",p_value->Y[0][1][r][c]);
        }printf("\n");
    }*/
}

int IDCT(YCbCr *p_value, my2DDCT *p_dct) {
    // f[r,c] -> F[u,v]
    int i, j, u, v, r, c;
    int H = 8, W = 8;

    for (i = 0; i < p_value->count_H; i++) {
        for (j = 0; j < p_value->count_W; j++) {
            for (r = 0; r < p_dct->H; r++) {
                for (c = 0; c < p_dct->W; c++) {
                    p_value->Y[i][j][r][c] = 0.0;
                    p_value->Cb[i][j][r][c] = 0.0;
                    p_value->Cr[i][j][r][c] = 0.0;
                    for (u = 0; u < p_dct->H; u++) {
                        for (v = 0; v < p_dct->W; v++) {
                            // 计算 DCT 基函数
                            p_dct->basis[u][v][r][c] = beta(u) * beta(v) * cos(PI * u * (2.0 * r + 1.0) / 2.0 / p_dct->H) * cos(PI * v * (2.0 * c + 1.0) / 2.0 / p_dct->W);
                            // 计算 Y 基函数的 IDCT 变换
                            p_value->Y[i][j][r][c] += (p_dct->Y_basis[i][j][u][v]) * p_dct->basis[u][v][r][c];
                            p_value->Cb[i][j][r][c] += (p_dct->Cb_basis[i][j][u][v]) * p_dct->basis[u][v][r][c];
                            p_value->Cr[i][j][r][c] += (p_dct->Cr_basis[i][j][u][v]) * p_dct->basis[u][v][r][c];
                            //if( i==1 && j==1) printf("%f\n",p_dct->Y_basis[i][j][r][c]);
                            
                        }
                    }

                    // 乘以逆变换的缩放因子
                    //if( i==1 && j==1) printf("before :%f\t",p_value->Y[i][j][r][c]);
                    p_value->Y[i][j][r][c] = p_value->Y[i][j][r][c] * sqrt(4.0 / H / W);
                    p_value->Cb[i][j][r][c] = p_value->Cb[i][j][r][c] * sqrt(4.0 / H / W);
                    p_value->Cr[i][j][r][c] = p_value->Cr[i][j][r][c] * sqrt(4.0 / H / W);
                    //if( i==1 && j==1) printf("%f\n",p_value->Y[i][j][r][c]);
                }
            }
        }
    }

    // 输出调试信息
    /*printf("decode DCT\n");
    for(r=0;r<8;r++){
        for(c=0;c<8;c++){
            printf("%f\t",p_value->Y[0][1][r][c]);
        }printf("\n");
    }*/
}

/* 以 ascii 儲存的channel  的 quantization table & 將值存到cut8_8_qt struct中*/
int quantization_table(YCbCr *p_value)
{
    int i,j;
    int H=8,W=8;
    int Y_quantization_table[8][8] = {
    {16,11,10,16,24,40,51,61},
    {12,12,14,19,26,58,60,55},
    {14,13,16,24,40,57,69,56},
    {14,17,22,29,51,87,80,62},
    {18,22,37,56,68,109,103,77},
    {24,35,55,64,81,104,113,92},
    {49,64,78,87,103,121,120,101},
    {72,92,95,98,112,100,103,99}
    };

    int CbCr_quantization_table[8][8]= {
    {17,18,24,47,99,99,99,99},
    {18,21,26,66,99,99,99,99},
    {24,26,56,99,99,99,99,99},
    {47,66,99,99,99,99,99,99},
    {99,99,99,99,99,99,99,99},
    {99,99,99,99,99,99,99,99},
    {99,99,99,99,99,99,99,99},
    {99,99,99,99,99,99,99,99}
    };

    for(i=0;i<H;i++){
        for(j=0;j<W;j++){
            p_value->Y_qt[i][j]= Y_quantization_table[i][j];
            p_value->Cb_qt[i][j]= CbCr_quantization_table[i][j];
            p_value->Cr_qt[i][j]= CbCr_quantization_table[i][j];
        }
    }
    return 1;
}

/* DCT後的值quantization */
int quantize(my2DDCT *p_dct, YCbCr *p_value){
    int i,j,u,v;
    int H=8,W=8;
   // printf("encode quantize\n");
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    p_dct->Y[i][j][u][v] = (int)round((p_dct->Y_basis[i][j][u][v])/((float)p_value->Y_qt[u][v]));
                    p_dct->Cb[i][j][u][v] = (int)round((p_dct->Cb_basis[i][j][u][v])/((float)p_value->Cb_qt[u][v]));
                    p_dct->Cr[i][j][u][v] = (int)round((p_dct->Cr_basis[i][j][u][v])/((float)p_value->Cr_qt[u][v]));
                    //calculate error
                    p_dct->eY[i][j][u][v] = (p_dct->Y_basis[i][j][u][v]) - ((p_value->Y_qt[u][v])*(p_dct->Y[i][j][u][v]));
                    p_dct->eCb[i][j][u][v] = (p_dct->Cb_basis[i][j][u][v]) - ((p_value->Cb_qt[u][v])*(p_dct->Cb[i][j][u][v]));
                    p_dct->eCr[i][j][u][v] = (p_dct->Cr_basis[i][j][u][v]) - ((p_value->Cr_qt[u][v])*(p_dct->Cr[i][j][u][v]));
                    //if(i==0 && j==0) printf("%f\t",p_dct->Cb_basis[0][0][u][v]); //DEBUG
                }
                //if(i==0 && j==0) printf("\n");
            }
        }
    }
    /*for(u=0;u<H;u++){
        for(v=0;v<W;v++){
            printf("%d\t",p_dct->Y[0][1][u][v]);
            }
        }*/
    return 1;
}

int I_quantize(YCbCr *p_value,my2DDCT *p_dct){
    int i,j,u,v ;
    int H=8,W=8;
    //printf("decode quantize\n");
    /*for(u=0;u<H;u++){
        for(v=0;v<W;v++){
            printf("%d\t",p_dct->Y[0][1][u][v]);
            }
        }*/
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    p_dct->Y_basis[i][j][u][v] = (float)p_dct->Y[i][j][u][v] *((float)p_value->Y_qt[u][v])+p_dct->eY[i][j][u][v];
                    p_dct->Cb_basis[i][j][u][v] = (float)p_dct->Cb[i][j][u][v] * ((float)p_value->Cb_qt[u][v])+p_dct->eCb[i][j][u][v];
                    p_dct->Cr_basis[i][j][u][v] = (float)p_dct->Cr[i][j][u][v] * ((float)p_value->Cr_qt[u][v])+p_dct->eCr[i][j][u][v];
                    //if(i==1 && j==1) printf("%f\t",p_dct->Y_basis[1][1][u][v]);
                }
                //if(i==0 && j==0) printf("\n");
            }
        }
    }
    return 1;
}
/* calculate SQNR*/
int SQNR(my2DDCT *p_dct,YCbCr *p_value){
    int i,j,u,v;
    float sqnr_Y=0.0,sqnr_eY = 0.0 ;
    float sqnr_Cb=0.0,sqnr_eCb = 0.0 ;
    float sqnr_Cr=0.0,sqnr_eCr = 0.0 ;
    float SQNR_Y=0.0,SQNR_Cb=0.0,SQNR_Cr=0.0;
    int H=8,W=8;
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    sqnr_Y = sqnr_Y + pow(p_dct->Y_basis[i][j][u][v],2) ;
                    sqnr_Cb = sqnr_Cb +pow(p_dct->Cb_basis[i][j][u][v],2) ;
                    sqnr_Cr = sqnr_Cr + pow(p_dct->Cr_basis[i][j][u][v],2) ;
                    sqnr_eY = sqnr_eY + pow(p_dct->eY[i][j][u][v],2) ;
                    sqnr_eCb = sqnr_eCb + pow(p_dct->eCb[i][j][u][v],2) ;
                    sqnr_eCr = sqnr_eCr + pow(p_dct->eCr[i][j][u][v],2) ;             
                }

            }
        }
    }
    SQNR_Y = 10 * log10(sqnr_Y /sqnr_eY );
    SQNR_Cb = 10 * log10(sqnr_Cb /sqnr_eCb );
    SQNR_Cr = 10 * log10(sqnr_Cr /sqnr_eCr );
    //printf("%f\t%f\n",sqnr_sum,sqnr_esum);
    printf("SQNR_Y: %f\n",SQNR_Y);
    printf("SQNR_Cb: %f\n",SQNR_Cb);
    printf("SQNR_Cr: %f\n",SQNR_Cr);
    return 1;
    
}

/* DPCM DC*/
int DPCM(my2DDCT *p_dct,YCbCr *p_value){
    int i,j;
    for(i=0;i<p_value->count_H;i++){
        for(j=p_value->count_W-1;j>=1;j--){
            p_dct->Y[i][j][0][0] = p_dct->Y[i][j][0][0]-p_dct->Y[i][j-1][0][0];
            p_dct->Cb[i][j][0][0] = p_dct->Cb[i][j][0][0]-p_dct->Cb[i][j-1][0][0];
            p_dct->Cr[i][j][0][0] = p_dct->Cr[i][j][0][0]-p_dct->Cr[i][j-1][0][0];
        }  
    }
    
    for(i=p_value->count_H-1;i>=1;i--){
        p_dct->Y[i][0][0][0] = p_dct->Y[i][0][0][0]-p_dct->Y[i-1][0][0][0];
        p_dct->Cb[i][0][0][0] = p_dct->Cb[i][0][0][0]-p_dct->Cb[i-1][0][0][0];
        p_dct->Cr[i][0][0][0] = p_dct->Cr[i][0][0][0]-p_dct->Cr[i-1][0][0][0];
        
    }
    return 1;
}

/* zigzag */
int zigzag(my2DDCT *p_dct,YCbCr *p_value){
    int i,j,u,v;
    int H=8,W=8;
    int zz_matrix[8][8] = {
    {0, 1, 5, 6, 14, 15, 27, 28},
    {2, 4, 7, 13, 16, 26, 29, 42},
    {3, 8, 12, 17, 25, 30, 41, 43},
    {9, 11, 18, 24, 31, 40, 44, 53},
    {10, 19, 23, 32, 39, 45, 52, 54},
    {20, 22, 33, 38, 46, 51, 55, 60},
    {21, 34, 37, 47, 50, 56, 59, 61},
    {35, 36, 48, 49, 57, 58, 62, 63}
    };
    
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    p_value->Y_zigzag[i][j][zz_matrix[u][v]]=p_dct->Y[i][j][u][v] ;
                    p_value->Cb_zigzag[i][j][zz_matrix[u][v]]=p_dct->Cb[i][j][u][v] ;
                    p_value->Cr_zigzag[i][j][zz_matrix[u][v]]=p_dct->Cr[i][j][u][v] ;
                }
                
            }
        }
    }
    /*printf("Before zigzag Y:\n");
    for(i=0;i<8;i++){
        for(j=0;j<8;j++){
            printf("%d\t",p_dct->Y[0][12][i][j]);
        }
        printf("\n");  
    }
    puts("\n");
    printf("After zigzag Y:\n");
    for(i=0;i<64;i++){
        printf("%d\t",p_value->Y_zigzag[0][12][i]);
        if((i+1)%8==0) printf("\n");
    }
    printf("\n");*/
    return 1;
}

/* Run Length Encoding (RLE) */
int RLE1(YCbCr *p_value){
    int i, j, k;
    int H = 8, W = 8;
    int nzero, nzero_count;
    
    for (i = 0; i < p_value->count_H; i++) {
        for (j = 0; j < p_value->count_W; j++) {
            nzero = 0;
            nzero_count = 0;
            k = 0;
            while (k < H * W) {
                //<nzero , zigzag[][][]>
                if (p_value->Y_zigzag[i][j][k] != 0) {
                    // Reallocate memory for RLE array
                    p_value->Y_rle[i][j] = (int*)realloc(p_value->Y_rle[i][j], (2 * nzero_count + 2) * sizeof(int));
                    p_value->Y_rle[i][j][2 * nzero_count] = nzero;
                    p_value->Y_rle[i][j][2 * nzero_count + 1] = p_value->Y_zigzag[i][j][k];
                    nzero = 0;
                    nzero_count++;//存放在哪個陣列位置
                } else if (p_value->Y_zigzag[i][j][k] == 0) {
                    nzero++;
                }
                k++;
            }
            //補0 0方便後續判斷
            p_value->Y_rle[i][j] = (int*)realloc(p_value->Y_rle[i][j], (2 * nzero_count + 2) * sizeof(int));
            p_value->Y_rle[i][j][2 * nzero_count] = 0;
            p_value->Y_rle[i][j][2 * nzero_count + 1] = 0;
        }
    }
    /*printf("Y_rle[0][12]:\n");
    for(k=0;k<14;k++){
        printf("%d\t",p_value->Y_rle[0][12][k]);
    }
    puts("\n");
    */  
    for (i = 0; i < p_value->count_H; i++) {
        for (j = 0; j < p_value->count_W; j++) {
            nzero = 0;
            nzero_count = 0;
            k = 0;
            while (k < H * W) {
                if (p_value->Cb_zigzag[i][j][k] != 0) {
                    // Reallocate memory for RLE array
                    p_value->Cb_rle[i][j] = (int*)realloc(p_value->Cb_rle[i][j], (2 * nzero_count + 2) * sizeof(int));
                    p_value->Cb_rle[i][j][2 * nzero_count] = nzero;
                    p_value->Cb_rle[i][j][2 * nzero_count + 1] = p_value->Cb_zigzag[i][j][k];
                    nzero = 0;
                    nzero_count++;
                } else if (p_value->Cb_zigzag[i][j][k] == 0) {
                    nzero++;
                }
                k++;
            }
            // Reallocate memory for the end of block
            p_value->Cb_rle[i][j] = (int*)realloc(p_value->Cb_rle[i][j], (2 * nzero_count + 2) * sizeof(int));
            p_value->Cb_rle[i][j][2 * nzero_count] = 0;
            p_value->Cb_rle[i][j][2 * nzero_count + 1] = 0;
        }
    }
    for (i = 0; i < p_value->count_H; i++) {
        for (j = 0; j < p_value->count_W; j++) {
            nzero = 0;
            nzero_count = 0;
            k = 0;
            while (k < H * W) {
                if (p_value->Cr_zigzag[i][j][k] != 0) {
                    // Reallocate memory for RLE array
                    p_value->Cr_rle[i][j] = (int*)realloc(p_value->Cr_rle[i][j], (2 * nzero_count + 2) * sizeof(int));
                    p_value->Cr_rle[i][j][2 * nzero_count] = nzero;
                    p_value->Cr_rle[i][j][2 * nzero_count + 1] = p_value->Cr_zigzag[i][j][k];
                    nzero = 0;
                    nzero_count++;
                } else if (p_value->Cr_zigzag[i][j][k] == 0) {
                    nzero++;
                }
                k++;
            }
            // Reallocate memory for the end of block
            p_value->Cr_rle[i][j] = (int*)realloc(p_value->Cr_rle[i][j], (2 * nzero_count + 2) * sizeof(int));
            p_value->Cr_rle[i][j][2 * nzero_count] = 0;
            p_value->Cr_rle[i][j][2 * nzero_count + 1] = 0;
        }
    }
    //demo2 之後
    size_t n ;
    size_t rle_code_size = 0;
    size_t block_size = 0 ;
    size_t one_block_size = 0;
    
    for ( i = 0; i <p_value->count_H; i++) {
        for (j = 0; j < p_value->count_W; j++) {
            // 取得合併陣列空間大小
            n =0;//skip
            while(((p_value->Y_rle[i][j][n])!= 0 ) || ((p_value->Y_rle[i][j][n+1])!=0)){
                n++;
            }
                n+=2;
                block_size += n ;
            }
        }
    //printf("Y size: %zu\n",block_size);
    p_value->Y_block_size = block_size ;
        
    // 僵值assign給該陣列
    p_value->Y_rle_code = (int*)realloc(p_value->Y_rle_code, block_size * sizeof(int));
   
    for ( i = 0; i < p_value->count_H; i++) {
        for ( j = 0; j < p_value->count_W; j++) {
                n=0;
            while(((p_value->Y_rle[i][j][n])!= 0 )|| ((p_value->Y_rle[i][j][n+1])!=0)){
                n++;
            }
                n+=2;
                one_block_size = n ;
            for ( k = 0; k < one_block_size; k++) {
                p_value->Y_rle_code[rle_code_size++] = p_value->Y_rle[i][j][k];
            }
        }
    }
    
    block_size = 0 ;
    rle_code_size = 0;
    for ( i = 0; i < p_value->count_H; i++) {
        for ( j = 0; j < p_value->count_W; j++) {
            // 获取当前块的大小
            n =0;
            while(p_value->Cb_rle[i][j][n]!=0 || p_value->Cb_rle[i][j][n+1]!=0){
                n++;
            }
                n+=2;
                block_size += n ;
            }
        }
        
    p_value->Cb_block_size = block_size ;

    p_value->Cb_rle_code = (int*)realloc(p_value->Cb_rle_code, block_size * sizeof(int));
   
    for ( i = 0; i < p_value->count_H; i++) {
        for ( j = 0; j < p_value->count_W; j++) {
                n=0;
            while(p_value->Cb_rle[i][j][n]!=0 || p_value->Cb_rle[i][j][n+1]!=0){
                n++;
            }
                n+=2;
                one_block_size = n ;
            for ( k = 0; k < one_block_size; k++) {
                p_value->Cb_rle_code[rle_code_size++] = p_value->Cb_rle[i][j][k];
            }
        }
    }
    //printf("CB size: %zu\n",rle_code_size);

    block_size = 0 ;
    rle_code_size = 0;
    for ( i = 0; i < p_value->count_H; i++) {
        for ( j = 0; j < p_value->count_W; j++) {
            
            n =0;
            while(p_value->Cr_rle[i][j][n]!=0 || p_value->Cr_rle[i][j][n+1]!=0){
                n++;
            }
                n+=2;
                block_size += n ;
            }
        }
        //printf("Cr size: %zu\n",block_size);
        p_value->Cr_block_size = block_size ;
        
    p_value->Cr_rle_code = (int*)realloc(p_value->Cr_rle_code, block_size * sizeof(int));
    rle_code_size = 0;
    for ( i = 0; i < p_value->count_H; i++) {
        for (j = 0; j < p_value->count_W; j++) {
                n=0;
            while(p_value->Cr_rle[i][j][n]!=0 || p_value->Cr_rle[i][j][n+1]!=0){
                n++;
            }
                n+=2;
                one_block_size = n ;
            for ( k = 0; k < one_block_size; k++) {
                p_value->Cr_rle_code[rle_code_size++] = p_value->Cr_rle[i][j][k];
            }
        }
    }
    return 1;
}
/* RLC 整併(0 3 4 5 3 4 -> 0:1,3:2,4:2,5:1)*/
int RLE2(YCbCr *p_value, RLE *p_rle){
    p_rle->Y_num =calloc(p_value->Y_block_size,sizeof(int)) ;
    p_rle->Y_count = calloc(p_value->Y_block_size,sizeof(int));
    p_rle->Cb_num =calloc(p_value->Cb_block_size,sizeof(int)) ;
    p_rle->Cb_count = calloc(p_value->Cb_block_size,sizeof(int));
    p_rle->Cr_num =calloc(p_value->Cr_block_size,sizeof(int)) ;
    p_rle->Cr_count = calloc(p_value->Cr_block_size,sizeof(int));
    int i,j ;
    int Y_count=1;
    int foundvalue1,foundvalue2;
    p_rle->Y_num[0] = p_value->Y_rle_code[0] ;
    p_rle->Y_count[0] = 1 ;
    for(i=1;i<p_value->Y_block_size;i++){
        foundvalue1=0,foundvalue2=0;
        for(j=0 ;j<Y_count;j++){
            if(p_value->Y_rle_code[i]!=p_rle->Y_num[j]){
                foundvalue1 = 1;
            }else if(p_value->Y_rle_code[i]==p_rle->Y_num[j]){
                foundvalue2 = 1;
                p_rle->Y_count[j]=p_rle->Y_count[j]+1;
            }
        }
        if(foundvalue1==1 &&foundvalue2==0){
            Y_count++;
            p_rle->Y_num =(int*)realloc(p_rle->Y_num, Y_count * sizeof(int));
            p_rle->Y_count = (int*)realloc(p_rle->Y_count, Y_count * sizeof(int));
            p_rle->Y_num[Y_count-1] = p_value->Y_rle_code[i] ;
            p_rle->Y_count[Y_count-1] = 1 ;
        }         
    }
    p_rle->Y_len = Y_count ;

    int Cb_count = 1;
    p_rle->Cb_num[0] = p_value->Cb_rle_code[0] ;
    p_rle->Cb_count[0] = 1 ;
    for(i=1;i<p_value->Cb_block_size;i++){
        foundvalue1=0,foundvalue2=0;
        for(j=0 ;j<Cb_count;j++){
            if(p_value->Cb_rle_code[i]!=p_rle->Cb_num[j]){
                foundvalue1 = 1;
            }else if(p_value->Cb_rle_code[i]==p_rle->Cb_num[j]){
                foundvalue2 = 1;
                p_rle->Cb_count[j]=p_rle->Cb_count[j]+1;
            }
        }
        if(foundvalue1==1 &&foundvalue2==0){
            Cb_count++;
            p_rle->Cb_num =(int*)realloc(p_rle->Cb_num, Cb_count * sizeof(int));
            p_rle->Cb_count = (int*)realloc(p_rle->Cb_count, Cb_count * sizeof(int));
            p_rle->Cb_num[Cb_count-1] = p_value->Cb_rle_code[i] ;
            p_rle->Cb_count[Cb_count-1] = 1 ;
        }         
    }
    p_rle->Cb_len = Cb_count ;
    int Cr_count = 1;
    p_rle->Cr_num[0] = p_value->Cr_rle_code[0] ;
    p_rle->Cr_count[0] = 1 ;
    for(i=1;i<p_value->Cr_block_size;i++){
        foundvalue1=0,foundvalue2=0;
        for(j=0 ;j<Cr_count;j++){
            if(p_value->Cr_rle_code[i]!=p_rle->Cr_num[j]){
                foundvalue1 = 1;
            }else if(p_value->Cr_rle_code[i]==p_rle->Cr_num[j]){
                foundvalue2 = 1;
                p_rle->Cr_count[j]=p_rle->Cr_count[j]+1;
            }
        }
        if(foundvalue1==1 &&foundvalue2==0){
            Cr_count++;
            p_rle->Cr_num =(int*)realloc(p_rle->Cr_num, Cr_count * sizeof(int));
            p_rle->Cr_count = (int*)realloc(p_rle->Cr_count, Cr_count * sizeof(int));
            p_rle->Cr_num[Cr_count-1] = p_value->Cr_rle_code[i] ;
            p_rle->Cr_count[Cr_count-1] = 1 ;
        }         
    }
    p_rle->Cr_len = Cr_count ;
    /* bubble sort */
    int tmp_count,tmp_num ;
    for(int i=0; i<Y_count-1; i++){
        for(int j=0; j<Y_count-1-i; j++){
            if(p_rle->Y_count[j] > p_rle->Y_count[j+1]){
                tmp_count = p_rle->Y_count[j]  ;
                p_rle->Y_count[j] = p_rle->Y_count[j+1];
                p_rle->Y_count[j+1] = tmp_count ;
                tmp_num = p_rle->Y_num[j]  ;
                p_rle->Y_num[j] = p_rle->Y_num[j+1];
                p_rle->Y_num[j+1] = tmp_num ;
            }	
        }
    }
    for(int i=0; i<Cb_count-1; i++){
        for(int j=0; j<Cb_count-1-i; j++){
            if(p_rle->Cb_count[j] > p_rle->Cb_count[j+1]){
                tmp_count = p_rle->Cb_count[j]  ;
                p_rle->Cb_count[j] = p_rle->Cb_count[j+1];
                p_rle->Cb_count[j+1] = tmp_count ;
                tmp_num = p_rle->Cb_num[j]  ;
                p_rle->Cb_num[j] = p_rle->Cb_num[j+1];
                p_rle->Cb_num[j+1] = tmp_num ;
            }	
        }
    }
    for( i=0; i<Cr_count-1; i++){
        for( j=0; j<Cr_count-1-i; j++){
            if(p_rle->Cr_count[j] > p_rle->Cr_count[j+1]){
                tmp_count = p_rle->Cr_count[j]  ;
                p_rle->Cr_count[j] = p_rle->Cr_count[j+1];
                p_rle->Cr_count[j+1] = tmp_count ;
                tmp_num = p_rle->Cr_num[j]  ;
                p_rle->Cr_num[j] = p_rle->Cr_num[j+1];
                p_rle->Cr_num[j+1] = tmp_num ;
            }	
        }
    }
    /*puts("\n");
    printf("After merge Y:\n");
    for(i=0;i<Y_count;i++){
        printf("%d\t%d\n",p_rle->Y_num[i],p_rle->Y_count[i]);
    }
    puts("\n");
    */
    return 1;	
    
}

//initial 

/* initialize 2D-DCT bases */
int my2DDCT_initial(my2DDCT *p_dct, YCbCr *p_value)
{
    size_t i, j, k;
    size_t H=8,W=8;
	p_dct->H = H;
    p_dct->W = W;
    int u,v,r,c;
    //p_dct->basis[H][W][H][W];
    p_dct->basis = (float ****) calloc(p_dct->H, sizeof(float ***));
    for(i=0;i<p_dct->H;i++) {
        p_dct->basis[i] = (float ***) calloc(p_dct->W, sizeof(float **));
        for(j=0;j<p_dct->W;j++) {
            p_dct->basis[i][j] = (float **) calloc(p_dct->H, sizeof(float *));
            for(k=0;k<p_dct->H;k++) {
                p_dct->basis[i][j][k] = (float *) calloc(p_dct->W, sizeof(float));
            }
        }
    }
    p_dct->Y_basis = (float ****) calloc(p_value->count_H, sizeof(float ***));
    for(i=0;i<p_value->count_H;i++) {
        p_dct->Y_basis[i] = (float ***) calloc(p_value->count_W, sizeof(float **));
        for(j=0;j<p_value->count_W;j++) {
            p_dct->Y_basis[i][j] = (float **) calloc(p_dct->H, sizeof(float *));
            for(k=0;k<p_dct->H;k++) {
                p_dct->Y_basis[i][j][k] = (float *) calloc(p_dct->W, sizeof(float));
            }
        }
    }
    p_dct->Cb_basis = (float ****) calloc(p_value->count_H, sizeof(float ***));
    for(i=0;i<p_value->count_H;i++) {
        p_dct->Cb_basis[i] = (float ***) calloc(p_value->count_W, sizeof(float **));
        for(j=0;j<p_value->count_W;j++) {
            p_dct->Cb_basis[i][j] = (float **) calloc(p_dct->H, sizeof(float *));
            for(k=0;k<p_dct->H;k++) {
                p_dct->Cb_basis[i][j][k] = (float *) calloc(p_dct->W, sizeof(float));
            }
        }
    }
    p_dct->Cr_basis = (float ****) calloc(p_value->count_H, sizeof(float ***));
    for(i=0;i<p_value->count_H;i++) {
        p_dct->Cr_basis[i] = (float ***) calloc(p_value->count_W, sizeof(float **));
        for(j=0;j<p_value->count_W;j++) {
            p_dct->Cr_basis[i][j] = (float **) calloc(p_dct->H, sizeof(float *));
            for(k=0;k<p_dct->H;k++) {
                p_dct->Cr_basis[i][j][k] = (float *) calloc(p_dct->W, sizeof(float));
            }
        }
    }

    p_dct->Y = (int ****) calloc(p_value->count_H, sizeof(int ***));
    for(i=0;i<p_value->count_H;i++) {
        p_dct->Y[i] = (int ***) calloc(p_value->count_W, sizeof(int **));
        for(j=0;j<p_value->count_W;j++) {
            p_dct->Y[i][j] = (int **) calloc(H, sizeof(int *));
            for(k=0;k<p_dct->H;k++) {
                p_dct->Y[i][j][k] = (int *) calloc(W, sizeof(int));
            }
        }
    }
    p_dct->Cb = (int ****) calloc(p_value->count_H, sizeof(int ***));
    for(i=0;i<p_value->count_H;i++) {
        p_dct->Cb[i] = (int ***) calloc(p_value->count_W, sizeof(int **));
        for(j=0;j<p_value->count_W;j++) {
            p_dct->Cb[i][j] = (int **) calloc(H, sizeof(int *));
            for(k=0;k<p_dct->H;k++) {
                p_dct->Cb[i][j][k] = (int *) calloc(W, sizeof(int));
            }
        }
    }
    p_dct->Cr = (int ****) calloc(p_value->count_H, sizeof(int ***));
    for(i=0;i<p_value->count_H;i++) {
        p_dct->Cr[i] = (int ***) calloc(p_value->count_W, sizeof(int **));
        for(j=0;j<p_value->count_W;j++) {
            p_dct->Cr[i][j] = (int **) calloc(H, sizeof(int *));
            for(k=0;k<p_dct->H;k++) {
                p_dct->Cr[i][j][k] = (int *) calloc(W, sizeof(int));
            }
        }
    }
    p_dct->eY = (float ****) calloc(p_value->count_H, sizeof(float ***));
    for(i=0;i<p_value->count_H;i++) {
        p_dct->eY [i] = (float ***) calloc(p_value->count_W, sizeof(float **));
        for(j=0;j<p_value->count_W;j++) {
            p_dct->eY [i][j] = (float **) calloc(H, sizeof(float *));
            for(k=0;k<p_dct->H;k++) {
                p_dct->eY [i][j][k] = (float *) calloc(W, sizeof(float));
            }
        }
    }
    p_dct->eCb = (float ****) calloc(p_value->count_H, sizeof(float ***));
    for(i=0;i<p_value->count_H;i++) {
        p_dct->eCb [i] = (float ***) calloc(p_value->count_W, sizeof(float **));
        for(j=0;j<p_value->count_W;j++) {
            p_dct->eCb [i][j] = (float **) calloc(H, sizeof(float *));
            for(k=0;k<p_dct->H;k++) {
                p_dct->eCb [i][j][k] = (float *) calloc(W, sizeof(float));
            }
        }
    }
    p_dct->eCr = (float ****) calloc(p_value->count_H, sizeof(float ***));
    for(i=0;i<p_value->count_H;i++) {
        p_dct->eCr[i] = (float ***) calloc(p_value->count_W, sizeof(float **));
        for(j=0;j<p_value->count_W;j++) {
            p_dct->eCr[i][j] = (float **) calloc(H, sizeof(float *));
            for(k=0;k<p_dct->H;k++) {
                p_dct->eCr[i][j][k] = (float *) calloc(W, sizeof(float));
            }
        }
    }

	return 1;
}

/* initialize cut8_8_qt */
int cut8_8_initial(BMP_file *p_bmp,YCbCr *p_value)
{ 
    int i,j,k,x,y;
    int H=8,W=8;
    p_value->Y = (float ****) calloc(p_value->count_H, sizeof(float ***));
    for(i=0;i<p_value->count_H;i++) {
        p_value->Y[i] = (float ***) calloc(p_value->count_W, sizeof(float **));
        for(j=0;j<p_value->count_W;j++) {
            p_value->Y[i][j] = (float **) calloc(H, sizeof(float *));
            for(k=0;k<8;k++) {
                p_value->Y[i][j][k] = (float *) calloc(W, sizeof(float));
            }
        }
    }
    p_value->Cb = (float ****) calloc(p_value->count_H, sizeof(float ***));
    for(i=0;i<p_value->count_H;i++) {
        p_value->Cb[i] = (float ***) calloc(p_value->count_W, sizeof(float **));
        for(j=0;j<p_value->count_W;j++) {
            p_value->Cb[i][j] = (float **) calloc(H, sizeof(float *));
            for(k=0;k<8;k++) {
                p_value->Cb[i][j][k] = (float *) calloc(W, sizeof(float));
            }
        }
    }
    p_value->Cr = (float ****) calloc(p_value->count_H, sizeof(float ***));
    for(i=0;i<p_value->count_H;i++) {
        p_value->Cr[i] = (float ***) calloc(p_value->count_W, sizeof(float **));
        for(j=0;j<p_value->count_W;j++) {
            p_value->Cr[i][j] = (float **) calloc(H, sizeof(int *));
            for(k=0;k<8;k++) {
                p_value->Cr[i][j][k] = (float *) calloc(W, sizeof(float));
            }
        }
    }

    p_value->Y_qt = (int **) calloc(H, sizeof(int *));
    for(i=0;i<H;i++){
        p_value->Y_qt[i] = (int*) calloc(W, sizeof(int));
    }
    p_value->Cb_qt = (int **) calloc(H, sizeof(int *));
    for(i=0;i<H;i++){
        p_value->Cb_qt[i] = (int*) calloc(W, sizeof(int));
    }
    p_value->Cr_qt = (int **) calloc(H, sizeof(int *));
    for(i=0;i<H;i++){
        p_value->Cr_qt[i] = (int*) calloc(W, sizeof(int));
    }

    p_value->Y_zigzag = (int ***) calloc(p_value->count_H, sizeof(int **));
    for(i=0;i<p_value->count_H;i++){
        p_value->Y_zigzag[i] = (int**) calloc(p_value->count_W, sizeof(int*));
        for(j=0;j<p_value->count_W;j++){
            p_value->Y_zigzag[i][j] = (int*) calloc(H*W,sizeof(int));
        }
    }
    p_value->Cb_zigzag = (int ***) calloc(p_value->count_H, sizeof(int **));
    for(i=0;i<p_value->count_H;i++){
        p_value->Cb_zigzag[i] = (int**) calloc(p_value->count_W, sizeof(int*));
        for(j=0;j<p_value->count_W;j++){
            p_value->Cb_zigzag[i][j] = (int*) calloc(H*W,sizeof(int));
        }
    }
    p_value->Cr_zigzag = (int ***) calloc(p_value->count_H, sizeof(int **));
    for(i=0;i<p_value->count_H;i++){
        p_value->Cr_zigzag[i] = (int**) calloc(p_value->count_W, sizeof(int*));
        for(j=0;j<p_value->count_W;j++){
            p_value->Cr_zigzag[i][j] = (int*) calloc(H*W,sizeof(int));
        }
    }

    p_value->Y_rle = (int ***) calloc(p_value->count_H, sizeof(int **));
    for(i=0;i<p_value->count_H;i++){
        p_value->Y_rle[i] = (int**) calloc(p_value->count_W, sizeof(int*));
        for(j=0;j<p_value->count_W;j++){
            p_value->Y_rle[i][j]=NULL;
        }
    }
    
    p_value->Cb_rle = (int ***) calloc(p_value->count_H, sizeof(int **));
    for(i=0;i<p_value->count_H;i++){
        p_value->Cb_rle[i] = (int**) calloc(p_value->count_W, sizeof(int*));
        for(j=0;j<p_value->count_W;j++){
            p_value->Cb_rle[i][j]=NULL;
        }
    }
    p_value->Cr_rle = (int ***) calloc(p_value->count_H, sizeof(int **));
    for(i=0;i<p_value->count_H;i++){
        p_value->Cr_rle[i] = (int**) calloc(p_value->count_W, sizeof(int*));
        for(j=0;j<p_value->count_W;j++){
            p_value->Cr_rle[i][j]=NULL;
        }
    }
    p_value->Y_rle_code = NULL;
    p_value->Cb_rle_code = NULL;
    p_value->Cr_rle_code = NULL;
    

    return 1;
}

/* free 2D-DCT bases */
int my2DDCT_free(my2DDCT *p_dct,YCbCr *p_value)
{
    size_t i, j, k;    
    for(i=0;i<p_dct->H;i++) {        
        for(j=0;j<p_dct->W;j++) {
            for(k=0;k<p_dct->H;k++) {
                free(p_dct->basis[i][j][k]);
            }
            free(p_dct->basis[i][j]);
        }
        free(p_dct->basis[i]);
    }
    free(p_dct->basis);

    for(i=0;i<p_value->count_H;i++) {        
        for(j=0;j<p_value->count_W;j++) {
            for(k=0;k<p_dct->H;k++) {
                free(p_dct->Y_basis[i][j][k]);
                free(p_dct->Cb_basis[i][j][k]);
                free(p_dct->Cr_basis[i][j][k]);
                free(p_dct->Y[i][j][k]);
                free(p_dct->Cb[i][j][k]);
                free(p_dct->Cr[i][j][k]);
                free(p_dct->eY[i][j][k]);
                free(p_dct->eCb[i][j][k]);
                free(p_dct->eCr[i][j][k]);
            }
            free(p_dct->Y_basis[i][j]);
            free(p_dct->Cb_basis[i][j]);
            free(p_dct->Cr_basis[i][j]);
            free(p_dct->Y[i][j]);
            free(p_dct->Cb[i][j]);
            free(p_dct->Cr[i][j]);
            free(p_dct->eY[i][j]);
            free(p_dct->eCb[i][j]);
            free(p_dct->eCr[i][j]);
        }
        free(p_dct->Y_basis[i]);
        free(p_dct->Cb_basis[i]);
        free(p_dct->Cr_basis[i]);
        free(p_dct->Y[i]);
        free(p_dct->Cb[i]);
        free(p_dct->Cr[i]);
        free(p_dct->eY[i]);
        free(p_dct->eCb[i]);
        free(p_dct->eCr[i]);
    }
    free(p_dct->Y_basis);
    free(p_dct->Cb_basis);
    free(p_dct->Cr_basis);
    free(p_dct->Y);
    free(p_dct->Cb);
    free(p_dct->Cr);
    free(p_dct->eY);
    free(p_dct->eCb);
    free(p_dct->eCr);

    return 1;
}

/* free cut8_8_qt */
int cut8_8_free(YCbCr *p_value)
{
    size_t i, j, k;    
    for(i=0;i<p_value->count_H;i++) {        
        for(j=0;j<p_value->count_W;j++) {
            for(k=0;k<8;k++) {
                free(p_value->Y[i][j][k]);
                free(p_value->Cb[i][j][k]);
                free(p_value->Cr[i][j][k]);
            }
            free(p_value->Y[i][j]);
            free(p_value->Cb[i][j]);
            free(p_value->Cr[i][j]);
            free(p_value->Y_zigzag[i][j]);
            free(p_value->Cb_zigzag[i][j]);
            free(p_value->Cr_zigzag[i][j]);

            free(p_value->Y_rle[i][j]);
            free(p_value->Cb_rle[i][j]);
            free(p_value->Cr_rle[i][j]);
            
        }
        free(p_value->Y[i]);
        free(p_value->Cb[i]);
        free(p_value->Cr[i]);
        free(p_value->Y_zigzag[i]);
        free(p_value->Cb_zigzag[i]);
        free(p_value->Cr_zigzag[i]);
        free(p_value->Y_rle[i]);
        free(p_value->Cb_rle[i]);
        free(p_value->Cr_rle[i]);
    }
    free(p_value->Y);
    free(p_value->Cb);
    free(p_value->Cr);
    free(p_value->Y_zigzag);
    free(p_value->Cb_zigzag);
    free(p_value->Cr_zigzag);
    free(p_value->Y_rle);
    free(p_value->Cb_rle);
    free(p_value->Cr_rle);

    free(p_value->Y_rle_code);
    free(p_value->Cb_rle_code);
    free(p_value->Cr_rle_code);

    for(i=0;i<8;i++) {
        free(p_value->Y_qt[i]);
        free(p_value->Cb_qt[i]);
        free(p_value->Cr_qt[i]);
    } 
    free(p_value->Y_qt);
    free(p_value->Cb_qt);
    free(p_value->Cr_qt);  
 

    return 1;
}

//decode 1
int ckeck_size(BMP_file *p_bmp,my2DDCT *p_dct, YCbCr *p_value,char*fn_dim){

    FILE *file = fopen(fn_dim, "r");
    if (fscanf(file, "%d\t%d", &p_bmp->H, &p_bmp->W) != 2) {
        printf("讀取錯誤。\n");
        return 1;
    }
    p_bmp->f_H = p_bmp->H+(8-p_bmp->H%8);
    p_bmp->f_W = p_bmp->W+(8-p_bmp->W%8);
    p_dct->H = 8 ;
    p_dct->W = 8;
    p_value->count_H = p_bmp->f_H/8;
    p_value->count_W = p_bmp->f_W/8;
    int file_size = (p_bmp->H *p_bmp->W)*24/8;
    return file_size;
}
int input(YCbCr *p_value,char *Qt_Y,char *Qt_Cb,char *Qt_Cr){
    int i,j,u,v;
    int H=8,W=8;
    FILE* fn_QY = fopen(Qt_Y,"r");
    FILE* fn_QCb = fopen(Qt_Cb,"r");
    FILE* fn_QCr = fopen(Qt_Cr,"r");
    for(i=0;i<H;i++){
        for(j=0;j<W;j++){
            fscanf(fn_QY,"%d",&p_value->Y_qt[i][j]);
            fscanf(fn_QCb,"%d",&p_value->Cb_qt[i][j]);
            fscanf(fn_QCr,"%d",&p_value->Cr_qt[i][j]);
            p_value->Y_qt[i][j]=(int)p_value->Y_qt[i][j];
            p_value->Cb_qt[i][j]=(int)p_value->Cb_qt[i][j];
            p_value->Cr_qt[i][j]=(int)p_value->Cr_qt[i][j];
            //printf("%02x %02x %02x\t", p_value->Y_qt[i][j], p_value->Cb_qt[i][j], p_value->Cr_qt[i][j]); // DEBUG
        }
    }
    
    fclose(fn_QY);
    fclose(fn_QCb);
    fclose(fn_QCr);
    return 1;

}
int input_quantize(YCbCr *p_value, my2DDCT *p_dct,char *qF_Y,char* qF_Cb,char* qF_Cr,char* eF_Y,char *eF_Cb,char *eF_Cr){
    int i,j,u,v;
    int H=8,W=8;
    FILE *fp_Y = fopen(qF_Y, "rb");
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    fread(&p_dct->Y[i][j][u][v], sizeof(short), 1, fp_Y);
                }
            }
        }
    }

    fclose(fp_Y);
    FILE *fp_Cb = fopen(qF_Cb, "rb");
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    fread(&p_dct->Cb[i][j][u][v], sizeof(short), 1, fp_Cb);
                }
            }
        }
    }
    fclose(fp_Cb);
    FILE *fp_Cr = fopen(qF_Cr, "rb");
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    fread(&p_dct->Cr[i][j][u][v], sizeof(short), 1, fp_Cr);
                }
            }
        }
    }
    fclose(fp_Cr);
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    p_dct->Y[i][j][u][v] = p_dct->Y[i][j][u][v]-128;
                    p_dct->Cb[i][j][u][v] = p_dct->Cb[i][j][u][v]-128;
                    p_dct->Cr[i][j][u][v] = p_dct->Cr[i][j][u][v]-128;
                    //if(i==0 && j==0) printf("%f\t",p_dct->Cb_basis[0][0][u][v]); //DEBUG
                }
                //if(i==0 && j==0) printf("\n");
            }
        }
    }
    FILE *ep_Y = fopen(eF_Y, "rb");
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    fread(&p_dct->eY[i][j][u][v], sizeof(float), 1, ep_Y);
                }
            }
        }
    }
    fclose(ep_Y);
    FILE *ep_Cb = fopen(eF_Cb, "rb");
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    fread(&p_dct->eCb[i][j][u][v], sizeof(float), 1, ep_Cb);
                }
            }
        }
    }
    fclose(ep_Cb);

    FILE *ep_Cr = fopen(eF_Cr, "rb");
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    fread(&p_dct->eCr[i][j][u][v], sizeof(float), 1, ep_Cr);
                }
            }
        }
    }
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    p_dct->eY[i][j][u][v] = p_dct->eY[i][j][u][v]-128.0;
                    p_dct->eCb[i][j][u][v] = p_dct->eCb[i][j][u][v]-128.0;
                    p_dct->eCr[i][j][u][v] = p_dct->eCr[i][j][u][v]-128.0;
                    //if(i==0 && j==0) printf("%f\t",p_dct->Cb_basis[0][0][u][v]); //DEBUG
                }
                //if(i==0 && j==0) printf("\n");
            }
        }
    }
    fclose(ep_Cr);
    return 1;
}




int RGB_to_file(BMP_file *p_bmp,char *bmp_fn ,int file_size){
    int x,y;
    char skip_buf[3] = { 0, 0, 0 };
    FILE* outfile = fopen(bmp_fn,"wb");
    int death = file_size*8/(p_bmp->H*p_bmp->W);
    int diff=(p_bmp->W *p_bmp->H)*3-file_size;
    int skip=diff/3/p_bmp->H;

    p_bmp->header.identifier[0] = 'B';
    p_bmp->header.identifier[1] = 'M';
    p_bmp->header.filesize = file_size + 54;  // size + header
    p_bmp->header.reserved = 0;
    p_bmp->header.reserved2 = 0;
    p_bmp->header.bitmap_dataoffset = 54;  //file header + info header
    p_bmp->header.bitmap_headersize = 40;  // info header
    p_bmp->header.width = p_bmp->W  ;
    p_bmp->header.height = p_bmp->H;
    p_bmp->header.planes = 1;
    p_bmp->header.bits_perpixel =death ;  
    p_bmp->header.compression = 0;  
    p_bmp->header.bitmap_datasize = file_size;
    p_bmp->header.hresolution = 2835;
    p_bmp->header.vresolution = 2835;
    p_bmp->header.usedcolors = 0;
    p_bmp->header.importantcolors = 0;
    p_bmp->header.palette = 0;

	fwrite(&p_bmp->header.identifier, sizeof(short), 1, outfile);
	fwrite(&p_bmp->header.filesize, sizeof(int), 1, outfile);
	fwrite(&p_bmp->header.reserved, sizeof(short), 1, outfile);
	fwrite(&p_bmp->header.reserved2, sizeof(short), 1, outfile);
	fwrite(&p_bmp->header.bitmap_dataoffset, sizeof(int), 1, outfile);
	fwrite(&p_bmp->header.bitmap_headersize, sizeof(int), 1, outfile);
	fwrite(&p_bmp->header.width, sizeof(int), 1, outfile);
	fwrite(&p_bmp->header.height, sizeof(int), 1, outfile);
	fwrite(&p_bmp->header.planes, sizeof(short), 1, outfile);
	fwrite(&p_bmp->header.bits_perpixel, sizeof(short), 1, outfile);
	fwrite(&p_bmp->header.compression, sizeof(int), 1, outfile);
	fwrite(&p_bmp->header.bitmap_datasize, sizeof(int), 1, outfile);
	fwrite(&p_bmp->header.hresolution, sizeof(int), 1, outfile);
	fwrite(&p_bmp->header.vresolution, sizeof(int), 1, outfile);
	fwrite(&p_bmp->header.usedcolors, sizeof(int), 1, outfile);
	fwrite(&p_bmp->header.importantcolors, sizeof(int), 1, outfile);
   
	for (x = 0; x<p_bmp->H; x++){
		for (y = 0; y<p_bmp->W; y++){
			fwrite(&p_bmp->B[x][y], sizeof(char), 1, outfile);
			fwrite(&p_bmp->G[x][y], sizeof(char), 1, outfile);
			fwrite(&p_bmp->R[x][y], sizeof(char), 1, outfile);
			//if (x<10 && y<10) printf("%02x %02x %02x\t", p_bmp->R[x][y], p_bmp->G[x][y], p_bmp->B[x][y]); // DEBUG
		}
		if (skip != 0) { fwrite(skip_buf, skip, 1, outfile); }
		//if (x<10) printf("\n"); // DEBUG
	}
    
    fclose(outfile);
    return 1;
}



//output
int output_quantize_table(char *Qt_Y,char *Qt_Cb,char *Qt_Cr,YCbCr *p_value){
    FILE* fp_Y = fopen(Qt_Y,"w");
    int i,j;
    int H=8,W=8;
    for(i=0;i<H;i++){
        for(j=0;j<W;j++){
            fprintf(fp_Y,"%d\t",p_value->Y_qt[i][j]);
        }
        fprintf(fp_Y,"\n");
    }
    fclose(fp_Y);
    FILE* fp_Cb = fopen(Qt_Cb,"w");
    for(i=0;i<H;i++){
        for(j=0;j<W;j++){
            fprintf(fp_Cb,"%d\t",p_value->Cb_qt[i][j]);
        }
        fprintf(fp_Cb,"\n");
    }
    fclose(fp_Cb);
    FILE* fp_Cr = fopen(Qt_Cr,"w");
    for(i=0;i<H;i++){
        for(j=0;j<W;j++){
            fprintf(fp_Cr,"%d\t",p_value->Cr_qt[i][j]);
        }
        fprintf(fp_Cr,"\n");
    }
    fclose(fp_Cr);
    return 1;
}
int output_quantize(YCbCr *p_value, my2DDCT *p_dct,char *qF_Y,char* qF_Cb,char* qF_Cr){
    int i,j,u,v;
    int H=8,W=8;
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    p_dct->Y[i][j][u][v] = p_dct->Y[i][j][u][v]+128;
                    p_dct->Cb[i][j][u][v] = p_dct->Cb[i][j][u][v]+128;
                    p_dct->Cr[i][j][u][v] = p_dct->Cr[i][j][u][v]+128;
                    //if(i==0 && j==0) printf("%f\t",p_dct->Cb_basis[0][0][u][v]); //DEBUG
                }
                //if(i==0 && j==0) printf("\n");
            }
        }
    }

    FILE *fp_Y = fopen(qF_Y, "wb");
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    fwrite(&p_dct->Y[i][j][u][v], sizeof(short), 1, fp_Y);
                }
            }
        }
    }
    fclose(fp_Y);
    FILE *fp_Cb = fopen(qF_Cb, "wb");
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    fwrite(&p_dct->Cb[i][j][u][v], sizeof(short), 1, fp_Cb);
                }
            }
        }
    }
    fclose(fp_Cb);
    FILE *fp_Cr = fopen(qF_Cr, "wb");
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    fwrite(&p_dct->Cr[i][j][u][v], sizeof(short), 1, fp_Cr);
                }
            }
        }
    }
    fclose(fp_Cr);
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    p_dct->Y[i][j][u][v] = p_dct->Y[i][j][u][v]-128;
                    p_dct->Cb[i][j][u][v] = p_dct->Cb[i][j][u][v]-128;
                    p_dct->Cr[i][j][u][v] = p_dct->Cr[i][j][u][v]-128;
                    //if(i==0 && j==0) printf("%f\t",p_dct->Cb_basis[0][0][u][v]); //DEBUG
                }
                //if(i==0 && j==0) printf("\n");
            }
        }
    }
    return 1;
}

int output_quantize_error(YCbCr *p_value, my2DDCT *p_dct,char *eF_Y,char* eF_Cb,char* eF_Cr){
    int i,j,u,v;
    int H=8,W=8;
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    p_dct->eY[i][j][u][v] = p_dct->eY[i][j][u][v]+128.0;
                    p_dct->eCb[i][j][u][v] = p_dct->eCb[i][j][u][v]+128.0;
                    p_dct->eCr[i][j][u][v] = p_dct->eCr[i][j][u][v]+128.0;
                    //if(i==0 && j==0) printf("%f\t",p_dct->Cb_basis[0][0][u][v]); //DEBUG
                }
                //if(i==0 && j==0) printf("\n");
            }
        }
    }
    FILE *ep_Y = fopen(eF_Y, "wb");
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    fwrite(&p_dct->eY[i][j][u][v], sizeof(float), 1, ep_Y);
                }
            }
        }
    }
    fclose(ep_Y);
    FILE *ep_Cb = fopen(eF_Cb, "wb");
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    fwrite(&p_dct->eCb[i][j][u][v], sizeof(float), 1, ep_Cb);
                }
            }
        }
    }
    fclose(ep_Cb);
    FILE *ep_Cr = fopen(eF_Cr, "wb");
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    fwrite(&p_dct->eCr[i][j][u][v], sizeof(float), 1, ep_Cr);
                }
            }
        }
    }
    fclose(ep_Cr);
    for(i=0;i<p_value->count_H;i++){
        for(j=0;j<p_value->count_W;j++){
            for(u=0;u<H;u++){
                for(v=0;v<W;v++){
                    p_dct->eY[i][j][u][v] = p_dct->eY[i][j][u][v]-128;
                    p_dct->eCb[i][j][u][v] = p_dct->eCb[i][j][u][v]-128;
                    p_dct->eCr[i][j][u][v] = p_dct->eCr[i][j][u][v]-128;
                    //if(i==0 && j==0) printf("%f\t",p_dct->Cb_basis[0][0][u][v]); //DEBUG
                }
                //if(i==0 && j==0) printf("\n");
            }
        }
    }
    return 1;
}

int output_rle(BMP_file *p_bmp,YCbCr *p_value,char *rle_code){
    int i,j,k,m_Y=0,m_Cb=0,m_Cr=0,p,q;
    for ( i = 0; i < p_value->count_H; i++) {
        for (j = 0; j < p_value->count_W; j++) {
            k=0;
            while(p_value->Y_rle[i][j][k]!=0 || p_value->Y_rle[i][j][k+1]!=0){
                k+=2;
            }
            m_Y+=k;
            k=0;
            while(p_value->Cb_rle[i][j][k]!=0 || p_value->Cb_rle[i][j][k+1]!=0){
               
                k+=2;
            }
            m_Cb+=k;
            k=0;
      
            while(p_value->Cr_rle[i][j][k]!=0 || p_value->Cr_rle[i][j][k+1]!=0){
              
                k+=2;
            }
            m_Cr+=k;
        }
    }
    FILE* fn_rle = fopen(rle_code,"w");
    fprintf(fn_rle,"dim: %d\n",m_Cb+m_Cr+m_Y);
    for ( i = 0; i < p_value->count_H; i++) {
        for (j = 0; j < p_value->count_W; j++) {
            fprintf(fn_rle,"(%d , %d , Y ):",i,j);
            k=0;
            while(p_value->Y_rle[i][j][k]!=0 || p_value->Y_rle[i][j][k+1]!=0){
                fprintf(fn_rle,"<%d,%d> ",p_value->Y_rle[i][j][k],p_value->Y_rle[i][j][k+1]);
                k+=2;
            }
            fprintf(fn_rle,"\n");
            k=0;
            fprintf(fn_rle,"(%d , %d , Cb):",i,j);
            while(p_value->Cb_rle[i][j][k]!=0 || p_value->Cb_rle[i][j][k+1]!=0){
                fprintf(fn_rle,"<%d,%d> ",p_value->Cb_rle[i][j][k],p_value->Cb_rle[i][j][k+1]);
                k+=2;
            }
            fprintf(fn_rle,"\n");
            k=0;
            fprintf(fn_rle,"(%d , %d , Cr):",i,j);
            while(p_value->Cr_rle[i][j][k]!=0 || p_value->Cr_rle[i][j][k+1]!=0){
                fprintf(fn_rle,"<%d,%d> ",p_value->Cr_rle[i][j][k],p_value->Cr_rle[i][j][k+1]);
                k+=2;
            }
            fprintf(fn_rle,"\n");
        }
    }
    fclose(fn_rle);
}

