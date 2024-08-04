#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "myBMPLib.h"

/*read header*/
void _readheader(FILE* fp, BMP_header *x) {
	fread(&x->identifier, sizeof(x->identifier), 1, fp);
	fread(&x->filesize, sizeof(x->filesize), 1, fp);
	fread(&x->reserved, sizeof(x->reserved), 1, fp);
	fread(&x->reserved2, sizeof(x->reserved2), 1, fp);
	fread(&x->bitmap_dataoffset, sizeof(x->bitmap_dataoffset), 1, fp);
	fread(&x->bitmap_headersize, sizeof(x->bitmap_headersize), 1, fp);
	fread(&x->width, sizeof(x->width), 1, fp);
	fread(&x->height, sizeof(x->height), 1, fp);
	fread(&x->planes, sizeof(x->planes), 1, fp);
	fread(&x->bits_perpixel, sizeof(x->bits_perpixel), 1, fp);
	fread(&x->compression, sizeof(x->compression), 1, fp);
	fread(&x->bitmap_datasize, sizeof(x->bitmap_datasize), 1, fp);
	fread(&x->hresolution, sizeof(x->hresolution), 1, fp);
	fread(&x->vresolution, sizeof(x->vresolution), 1, fp);
	fread(&x->usedcolors, sizeof(x->usedcolors), 1, fp);
	fread(&x->importantcolors, sizeof(x->importantcolors), 1, fp);
}

/* load BMP given with file name (fn) */
int BMP_file_load_fn(char *fn, BMP_file *p_bmp)
{
    char tmp_data = 0;
    size_t i = 0, j = 0;
    int skip = 0;
    char skip_buf[3] = { 0, 0, 0 };
    FILE *fp = fopen(fn, "rb");// open file
    /* check if the file exist?*/
    if ( !fp ) {
        fprintf(stderr, "cannot open %s\n", fn);
        return 0;
    }

    strcpy(p_bmp->_filename, fn);
    
    // read header
    _readheader(fp, &(p_bmp->header));
    
    // data
    p_bmp->W = p_bmp->header.width;
    p_bmp->H = p_bmp->header.height;
    
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
    skip = (4 - (p_bmp->W * 3) % 4); //BMP字節為4，餘數不為4要補
    
	if (skip == 4) skip = 0;

    for (i = 0; i<p_bmp->H; i++){
		for (j = 0; j<p_bmp->W; j++){//BMP標準格式順序為B->G->R
			fread(&(p_bmp->B[i][j]), sizeof(unsigned char), 1, fp);
            fread(&(p_bmp->G[i][j]), sizeof(unsigned char), 1, fp);
            fread(&(p_bmp->R[i][j]), sizeof(unsigned char), 1, fp);
		}
		if (skip != 0) fread(skip_buf, skip, 1, fp);//4個字節最多可能有3個skip要跳過
	}
    /*W=3 3*3=9%4=1 , 4-1=3 
    B0 G0 R0 | B1 G1 R1 | B2 G2 R2 | value1 value2 */
    fclose(fp);
    return 1;
}



/* save BMP as text file */
int BMP_save_text(BMP_file *p_bmp, char *fn_R, char *fn_G, char *fn_B, char *fn_dim)
{
    FILE *fp_R = fopen(fn_R, "w");
    if (fp_R == NULL) {
        fprintf(stderr, "Don't open file %s\n", fp_R);
        return 1;
    }
    int i, j;
    for( i=0;i<p_bmp->H; i++) {
        for( j=0;j<p_bmp->W;j++) {
            fprintf(fp_R, "%d\t", p_bmp->R[i][j]);
        }
        fprintf(fp_R, "\n");
    }

    fclose(fp_R);

    FILE *fp_G = fopen(fn_G, "w");
    for( i=0;i<p_bmp->H; i++) {
        for( j=0;j<p_bmp->W;j++) {
            fprintf(fp_G, "%d\t", p_bmp->G[i][j]);
        }
        fprintf(fp_G, "\n");
    }

    fclose(fp_G);

    FILE *fp_B = fopen(fn_B, "w");
    for( i=0;i<p_bmp->H; i++) {
        for( j=0;j<p_bmp->W;j++) {
            fprintf(fp_B, "%d\t", p_bmp->B[i][j]);
        }
        fprintf(fp_B, "\n");
    }
    /*for (i = 0; i < 3 ;i++) {
        for (j = 0; j < 3; j++) {
                printf("%02x %02x %02x\t", p_bmp->R[i][j], p_bmp->G[i][j], p_bmp->B[i][j]); // DEBUG
            }
        printf("\n"); 
        }
    */

    fclose(fp_B);
    FILE *fp_size = fopen(fn_dim, "w");
    //8bits是否修正
    fprintf(fp_size,"%d",p_bmp->header.width * p_bmp->header.height*p_bmp->header.bits_perpixel/8);
    fclose(fp_size);


}

int output_size(BMP_file *p_bmp,char *fn_dim){
    FILE *fp_size = fopen(fn_dim, "w");
    fprintf(fp_size,"%d\t%d\t",p_bmp->header.height , p_bmp->header.width);
    fclose(fp_size);
    return 1;
}

int text_save_BMP(BMP_file *p_bmp,FILE* outfile, char *fn_R, char *fn_G, char *fn_B, int file_size){
	char skip_buf[3] = { 0, 0, 0 };
	int x, y,i,j;
    FILE *f_R = fopen(fn_R, "rb");
    
    char ch;
    p_bmp->H=0;
    p_bmp->W=0;
    // 逐字符讀取文件，計算行數和最大列數
    while ((ch = fgetc(f_R)) != EOF) {
        if (ch == '\n') 
            p_bmp->H=(p_bmp->H)+1;
        if(ch=='\t' && p_bmp->H==0)
            p_bmp->W=(p_bmp->W)+1; 
    }
    fclose(f_R);
    
    int death = file_size*8/(p_bmp->H*p_bmp->W);
    int diff=(p_bmp->W *p_bmp->H)*3-file_size;
    int skip=diff/3/p_bmp->H;
    //printf("H: %d , W: %d, D: %d, skip: %d\n",p_bmp->H,p_bmp->W,death,skip );//check H and W
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

    FILE *fp_R = fopen(fn_R,"rb");
    FILE *fp_G = fopen(fn_G,"rb");
    FILE *fp_B = fopen(fn_B,"rb");
    
    /* read text */
    for (i = 0; i < p_bmp->H; i++) {
        for (j = 0; j < p_bmp->W; j++) {
            fscanf(fp_B, "%d", &p_bmp->B[i][j]);
            fscanf(fp_G, "%d", &p_bmp->G[i][j]);
            fscanf(fp_R, "%d", &p_bmp->R[i][j]);
            p_bmp->B[i][j]=(int)p_bmp->B[i][j];
            p_bmp->G[i][j]=(int)p_bmp->G[i][j];
            p_bmp->R[i][j]=(int)p_bmp->R[i][j];
            //if (i<3 && j<3) printf("%02x %02x %02x\t", p_bmp->R[i][j], p_bmp->G[i][j], p_bmp->B[i][j]); // DEBUG
            }
            if (skip != 0) { fread(skip_buf, skip, 1, outfile); }
        //if (i<3) printf("\n");
        }
    
    fclose(fp_B);
    fclose(fp_G);
    fclose(fp_R);

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
			//if (x<3 && y<3) printf("%02x %02x %02x\t", p_bmp->R[x][y], p_bmp->G[x][y], p_bmp->B[x][y]); // DEBUG
		}
		if (skip != 0) { fwrite(skip_buf, skip, 1, outfile); }
		//if (x<3) printf("\n"); // DEBUG
	}

    
    
    
    fclose(outfile);
    return 1;
}

