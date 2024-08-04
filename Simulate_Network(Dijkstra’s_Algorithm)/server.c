#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/socket.h> 
#include <netinet/in.h>

int main (int argc, char **argv){
	struct sockaddr_in server, client;
	int sockfd , connfd;
	int clen;
	char buf[20];
	
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
        bzero(buf, sizeof(buf));
        read(connfd, buf, sizeof(buf));
        write(connfd, buf, strlen(buf));
        
        return 0;
}
