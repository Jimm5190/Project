INCLUDE = -I.

all: build_encoder build_decoder

build_encoder:
	gcc -O3 -o encoder.exe encoder.c myBMPLib.c my2DDCTLib.c -lm

build_decoder:
	gcc -O3 -o decoder.exe decoder.c myBMPLib.c my2DDCTLib.c -lm

demo0:
	./encoder.exe 0 KimberlyNCat.bmp R.txt G.txt B.txt dim.txt
	./decoder.exe 0 ResKimberlyNCat.bmp R.txt G.txt B.txt dim.txt
	diff KimberlyNCat.bmp ResKimberlyNCat.bmp
demo1:
	./encoder.exe 1 KimberlyNCat.bmp Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw qF_Cb.raw qF_Cr.raw eF_Y.raw eF_Cb.raw eF_Cr.raw
	./decoder.exe 1 QResKimberlyNCat.bmp Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw qF_Cb.raw qF_Cr.raw eF_Y.raw eF_Cb.raw eF_Cr.raw

demo2: 
	./encoder.exe 2 KimberlyNCat.bmp ascii rle_code.txt
demo3:
	./encoder.exe 3 KimberlyNCat.bmp Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw qF_Cb.raw qF_Cr.raw eF_Y.raw eF_Cb.raw eF_Cr.raw text.bmp
clean:
	rm -f encoder.exe decoder.exe R.txt G.txt B.txt dim.txt Qt_Y.txt Qt_Cb.txt Qt_Cr.txt qF_Y.raw qF_Cb.raw qF_Cr.raw eF_Y.raw eF_Cb.raw eF_Cr.raw rle_code.txt

.PHONY: build_encoder build_decoder demo0 clean
