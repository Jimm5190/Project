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
	packet_arrives[i]=gen_exp((double)lambda); //將lambda所形成變數依序儲存到陣列arrives作為第i個packet arrives 
	packet_arrives[0]=0; //宣告第一個packet arrives在第0秒時發生 
	//printf(" packet_arrives= %1f\n",packet_arrives[i]);    
}
for(i=0;i<value;i++){
    packet_leaves[i]=gen_exp((double)1);  //將1所形成變數依序儲存到陣列leaves作為第i個packet leaves
    //printf(" packet_leaves= %1f\n",packet_leaves[i]);
}
//printf(" packet 0 arrives= %1f\n", packet_arrives[0]);
for(i=1;i<value;i++){
	packet_arrives[i]=packet_arrives[i]+packet_arrives[i-1];//累加packet_arrives 
	//printf(" packet %d arrives= %1f\n", i ,packet_arrives[i]);
}
	
//printf(" packet 0 leaves= %1f\n", packet_leaves[0]);	
for(j=1;j<value;j++){ //累加packet_leaves 
	i=0;
	if(packet_arrives[j]>packet_leaves[i]){ //當累加的arrives值大於先前的leaves值時，第n個leaves值等於第n個arrives值加上原第n個leaves的值 
		packet_leaves[j]=packet_arrives[j]+packet_leaves[j];
		i=j; //宣告i=j避免比較完第n個值後，第n+1個值又跟第n個值去比 
		//printf(" packet %d leaves= %1f\n", i ,packet_leaves[i]);
	}
	else{ //當累加的arrives值小於先前的leaves值時，第n個leaves值等於第n個leaves值加上原第n-1個leaves的值 
		packet_leaves[j]=packet_leaves[j]+packet_leaves[j-1];
		i=j;
	    //printf(" packet %d leaves= %1f\n", i ,packet_leaves[i]);
	}
}

for(i=0;i<value;i++){ 
    	if(packet_arrives[i+1]<packet_leaves[i]){//如果arrives值小於leaves值時，num加一代表被discard數量 
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

srand((unsigned)time(NULL));//設定每次都隨機 

printf("Input lambda: ");
scanf("%lf", &lambda);
sim(lambda);

return 0;
}
