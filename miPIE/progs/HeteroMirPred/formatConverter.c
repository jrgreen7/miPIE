#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

int isNucleotide(char ch)
{
    return ch == 'A' || ch == 'C' || ch == 'G' || ch == 'U' || ch == 'N';
}

int main(int argc, char ** argv)
{
    int maxv = 0;
    if (argc != 4)
    {
	printf("Incorect number of parameters!\n");
	return 0;
    }
    
    FILE *f1, *f2;
    
    f1 = fopen(argv[1], "r");
    f2 = fopen(argv[2], "w");
    
    if (!f1)
    {
	printf("Cannot open file: %s\n", argv[1]);
	return 0;
    }
    
    if (!f2)
    {
	printf("Cannot open file: %s\n", argv[2]);
	return 0;
    }    

    int i;
    char s1[1024], s2[1024], s3[1024];
    fgets(s1, 1000, f1);
    while (1)
    {
	fgets(s2, 1000, f1);	
	memset(s3, 0, sizeof(s3));
	fgets(s3, 1000, f1);		
	int N;
	while (s3[0] && s3[0] != '>')
	{
            N = strlen(s2);
            while (!isNucleotide(s2[N - 1]))
		N--;	
	    if (isNucleotide(s3[0]))
	    {
		memcpy(&s2[N], s3, (strlen(s3) + 1) * sizeof(char));
		N = strlen(s2);
		while (!isNucleotide(s2[N - 1]))
		    N--;
		s2[N] = 0;
	    }
	    memset(s3, 0, sizeof(s3));
	    fgets(s3, 1000, f1);
	}		

        N = strlen(s2);
        while (!isNucleotide(s2[N - 1]))
	    N--;
	s2[N] = 0;
	
	maxv = MAX(maxv, strlen(s2));
	
	fprintf(f2, "%s %s\n", argv[3], s2);
	
	if (s3[0] != '>')
	    break;
    }
    
    printf("Maximum sequence length: %d\n", maxv);
    
    return 0;    
}