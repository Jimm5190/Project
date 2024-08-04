#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/socket.h> 
#include <netinet/in.h>
#include <stdlib.h> 


int path(int cost[10][10], int s, char obuf[100])
{
    int i,j,k,l,m,n,o,a ;
    int min[10];
    int node[10];
    int know[10];
    int answer[10];
    int order[10];
    int distance[10];
    
    for(i=0 ;i<10;i++){
    	distance[i]=10000;
    	min[i]=9999;
    	node[i] = -1 ;
    	order[i] = 10000;
    	answer[i] = 10000;
	}

    know[s]=1 ;
    node[0]=s ;
    for(i=0 ;i<6 ;i++)
    {
    	if(cost[s][i]<10000){
    		distance[i] = cost[s][i] ;
			node[i] = s ;}
        if(distance[i]<min[0]){
    	    min[0] = cost[s][i] ;
    	    distance[i]=min[0] ;
    	    k = i ;

			}		
	} 
	know[k]=1 ;
	
	for(i=0 ;i<6 ;i++){
		if(know[i]!=1 & cost[k][i]<10000){
			if((distance[k]+cost[k][i] ) < distance[i]){
    			        distance[i] = distance[k] +cost[k][i] ;
						node[i] = k ;
					}
		        }
		    } 
	   			
	for(i=0 ;i<6 ;i++){
			if(know[i]!=1){
    			if(distance[i]<min[1]){
    					min[1]=distance[i] ;
    					l = i ;
				}
			}
		}
	know[l] =1 ;

	for(i=0 ;i<6 ;i++){
		if(know[i]!=1 & cost[l][i]<10000){
			if((distance[l]+cost[l][i] ) < distance[i]){
    			        distance[i] = distance[l] +cost[l][i] ;
						node[i] = l ;}
		}
		}  			
	for(i=0 ;i<6 ;i++){
			if(know[i]!=1){
    			if(distance[i]<min[2]){
    					min[2]=distance[i] ;
    					m = i ;
				}
			}
		}
	know[m] =1 ;

	for(i=0 ;i<6 ;i++){
		if(know[i]!=1 & cost[m][i]<10000){
			if((distance[l]+cost[m][i] ) < distance[i]){
    			        distance[i] = distance[m] +cost[m][i] ;
						node[i] = m ;
						}
		            }
	        	}  			
	for(i=0 ;i<6 ;i++){
			if(know[i]!=1){
    			if(distance[i]<min[3]){
    					min[3]=distance[i] ;
    					n = i ;
				}
			}
		}
	know[n] =1 ;
    for(i=0 ;i<6 ;i++){
		if(know[i]!=1 & cost[n][i]<10000){
			if((distance[l]+cost[n][i] ) < distance[i]){
    			        distance[i] = distance[n] +cost[n][i] ;
						node[i] = n ;
					}
		        }
		    }  			
	for(i=0 ;i<6 ;i++){
			if(know[i]!=1){
    			if(distance[i]<min[4]){
    					min[4]=distance[i] ;
    					o = i ;
				}
			}
		}
	know[o] =1 ;
	
	for(i=0 ; i<6 ;i++){
	    if((distance[i]<cost[s][i])&&(cost[s][i]!=10000)){
	    	order[0] = s ;
	    	order[1] = i ;
	    	i = node[i] ;
	    	order[2] = i ;
	    	i = node[i] ;
	    	    if(node[i]==-1){
	    	        order[3] = s ;	
			        break ;
			}
			order[3] = i ;
	    	i = node[i] ;
			order[4] = s ;
	    	break;
		}
    }

    
    for(i = 0 ;i<4 ;i++){
    	answer[3-i] = order[i] ; 	
    	}
	
	if(order[4]<=5){
		for(i = 0 ;i<5 ;i++){
    	answer[4-i] = order[i] ;   	
	    }		
	}
    if(answer[4]>5){
        for (i = 0; i < 4 ; i++) {
            if (i<3) {
                sprintf(obuf + 2 * i, "%d,", answer[i]);          
        } 
		    else {
                sprintf(obuf + 2 * i, "%d", answer[i]);
            }
        }
	}
    if(answer[4]<=5){
        for (i = 0; i < 5 ; i++) {
            if (i<4) {
                sprintf(obuf + 2 * i, "%d,", answer[i]);
            } 
			else {
                sprintf(obuf + 2 * i, "%d", answer[i]);
            }
        }
    }
      if(answer[0] == answer[1]){
	    sprintf(obuf, "%s", "n");
      }
      

    return 0;  
}


int main(int argc, char **argv)
{
   
    int cost[10][10];
    int n,i,j,s,d;
    char obuf[100];
 
    for(i=0;i<10;i++)
    {
        for(j=0;j<10;j++)
        { 
           cost[i][j]=10000;
        }  
    }


    cost[0][1]=1;
    cost[0][3]=4;
    cost[1][0]=1;
    cost[1][2]=3;
    cost[1][4]=1;
    cost[2][1]=3;
    cost[2][4]=1;
    cost[2][5]=2;
    cost[3][0]=4;
    cost[3][4]=1;
    cost[4][1]=1;
    cost[4][2]=1;
    cost[4][3]=1;
    cost[5][2]=2;
    
    struct sockaddr_in server, client;
	int sockfd , connfd;
	int clen;
	char ibuf[20];
       

    server.sin_family = AF_INET;
	server.sin_addr.s_addr = INADDR_ANY;
	server.sin_port = htons(8888) ;
    
    sockfd = socket(AF_INET, SOCK_STREAM , 0);
	if(sockfd == -1){
		printf("Error: socket\n");
		return 0;
	}
	
	if(bind(sockfd,(struct sockaddr *)&server, sizeof(server))== -1)
        {
		printf("Error : bind\n");
		return 0 ;
	}
	
	if(listen(sockfd, 1)==-1)
	{
		printf("Error listen\n");
		return 0;
	}
	
	clen = sizeof(client);
	connfd = accept(sockfd , (struct sockaddr *)&client , &clen);
	if (connfd == -1){
		printf("Error: accept\n");
		return 0 ;
	}
	
	
        read(connfd, ibuf, sizeof(ibuf));
        s = atoi(ibuf);
        
    path(cost,s,obuf);
		
    write(connfd, obuf,strlen(obuf));

    return 0;
      
}
