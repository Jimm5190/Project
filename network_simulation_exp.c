#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define value 10000

double gen_exp(double lambda){
	double y, x;
	y=((double)rand()/((double)RAND_MAX+1.0));
	return x=-log(1-y)/lambda;
}


int sim(double lambda){
double packet_arrives[value];
double packet_leaves[value];
int i,j,num=0;
double F;

for(i=0;i<value;i++){
	packet_arrives[i]=gen_exp((double)lambda); //�Nlambda�ҧΦ��ܼƨ̧��x�s��}�Carrives�@����i��packet arrives 
	packet_arrives[0]=0; //�ŧi�Ĥ@��packet arrives�b��0��ɵo�� 
	//printf(" packet_arrives= %1f\n",packet_arrives[i]);    
}
for(i=0;i<value;i++){
    packet_leaves[i]=gen_exp((double)1);  //�N1�ҧΦ��ܼƨ̧��x�s��}�Cleaves�@����i��packet leaves
    //printf(" packet_leaves= %1f\n",packet_leaves[i]);
}
//printf(" packet 0 arrives= %1f\n", packet_arrives[0]);
for(i=1;i<value;i++){
	packet_arrives[i]=packet_arrives[i]+packet_arrives[i-1];//�֥[packet_arrives 
	//printf(" packet %d arrives= %1f\n", i ,packet_arrives[i]);
}
	
//printf(" packet 0 leaves= %1f\n", packet_leaves[0]);	
for(j=1;j<value;j++){ //�֥[packet_leaves 
	i=0;
	if(packet_arrives[j]>packet_leaves[i]){ //��֥[��arrives�Ȥj����e��leaves�ȮɡA��n��leaves�ȵ����n��arrives�ȥ[�W���n��leaves���� 
		packet_leaves[j]=packet_arrives[j]+packet_leaves[j];
		i=j; //�ŧii=j�קK�������n�ӭȫ�A��n+1�ӭȤS���n�ӭȥh�� 
		//printf(" packet %d leaves= %1f\n", i ,packet_leaves[i]);
	}
	else{ //��֥[��arrives�Ȥp����e��leaves�ȮɡA��n��leaves�ȵ����n��leaves�ȥ[�W���n-1��leaves���� 
		packet_leaves[j]=packet_leaves[j]+packet_leaves[j-1];
		i=j;
	    //printf(" packet %d leaves= %1f\n", i ,packet_leaves[i]);
	}
}

for(i=0;i<value;i++){ 
    	if(packet_arrives[i+1]<packet_leaves[i]){//�p�Garrives�Ȥp��leaves�ȮɡAnum�[�@�N��Qdiscard�ƶq 
    		num++; 
			}
		}
printf("num = %d\n",num);
		
F = num/(double)10000;

printf("Output: %1f\n",F);
	
return 0;
}


int main(){
double lambda;

srand((unsigned)time(NULL));//�]�w�C�����H�� 

printf("Input lambda: ");
scanf("%lf", &lambda);
sim(lambda);

return 0;
}
