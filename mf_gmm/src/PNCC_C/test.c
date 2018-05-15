#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
    if (strcmp(argv[1]+strlen(argv[1])-3,"wav")==0) 
        printf("Found wav\n");
    else if (strcmp(argv[1]+strlen(argv[1])-3,"sph")==0) 
        printf("Found sph\n");
    else printf("Nada\n");
    return 0; 
}
