#ifndef INCLUDE_myBMPLib
#define INCLUDE_myBMPLib

/*construct a structure of BMP header*/
typedef struct _BMP_header {
	char identifier[2]; // 0x0000
	unsigned int filesize; // 0x0002
	unsigned short reserved; // 0x0006
	unsigned short reserved2;
	unsigned int bitmap_dataoffset; // 0x000A
	unsigned int bitmap_headersize; // 0x000E
	unsigned int width; // 0x0012
	unsigned int height; // 0x0016
	unsigned short planes; // 0x001A
	unsigned short bits_perpixel; // 0x001C
	unsigned int compression; // 0x001E
	unsigned int bitmap_datasize; // 0x0022
	unsigned int hresolution; // 0x0026
	unsigned int vresolution; // 0x002A
	unsigned int usedcolors; // 0x002E
	unsigned int importantcolors; // 0x0032
	unsigned int palette; // 0x0036
} BMP_header;


/* structure for BMP file */
typedef struct _BMP_file{
	char _filename[1024];
    BMP_header header;
    int H;
	int W;
	int f_H;
	int f_W;
	unsigned char **B;// channel for Blue
	unsigned char **G;// channel for Green
	unsigned char **R;
	float **Y;
	float **Cb;
	float **Cr;
} BMP_file;

/* prototypes */

//demo0

/* load BMP given with file name (fn) */
int BMP_file_load_fn(char *fn, BMP_file *p_bmp);

/* save BMP as text file */
int BMP_save_text(BMP_file *p_bmp, char *fn_R, char *fn_G, char *fn_B, char *fn_dim);
/* return text to BMP */
int text_save_BMP(BMP_file *p_bmp, FILE* outfile, char *fn_R, char *fn_G, char *fn_B, int file_size);
/* output dim */
int output_size(BMP_file *p_bmp,char *fn_dim);
#endif