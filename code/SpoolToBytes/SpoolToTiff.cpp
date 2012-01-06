/*
 * Spool Header:  2*2752 Bytes = 5504   Bytes
 * Image Gap:     2*128  Bytes = 256    Bytes
 * Image Size:    2*512*512    = 524288 Bytes
 * Offset         10000 first file number
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#define HEADERSIZE 5504
#define IMAGESIZE 524288
#define GAPSIZE 256
#define TEXTSIZE 1024

int main(int argc, char **argv) 
{
  unsigned short image[IMAGESIZE];
  char gap[GAPSIZE];
  char textbuffer[TEXTSIZE];
  char *filebuffer;
  FILE *input_file, *output_file;
  char *spoolfilename = argv[1];
  int  MOD      = atoi( argv[2] );
  char A='a';

  if( argc != 3 )
  {
    fprintf(stderr, "Requires spool file name and modulus\n");
    fprintf(stderr, "SpoolToTiff /my/spool/file.spl 3\n");
    return 1;
  }
  
  printf("Splitting: %s\n", spoolfilename);
  
  input_file = fopen(spoolfilename,"rb");
  
  if (input_file == NULL)
  {
    fprintf(stderr, "Cannot open: %s\n", spoolfilename);
    return 2;
  }
  
  /* Read the entire spool into filebuffer */
  struct stat mystat;
  stat(spoolfilename,&mystat);
  filebuffer = (char*)malloc(mystat.st_size);
  size_t s = fread(filebuffer,1,mystat.st_size,input_file);
  fclose(input_file);
  if( s != mystat.st_size )
  {
    fprintf(stderr, "Size read (%lu) does not equal to file size (%lu)\n",
      s, mystat.st_size);
    return 3;
  }
  
  const int OFFSET = 10000;
  int numImages = (mystat.st_size - HEADERSIZE) / (GAPSIZE + IMAGESIZE);
  int count = -1;
  char *ptr = filebuffer + HEADERSIZE;
  printf( "Num Images: %d\n", numImages);
  for( int i = 0; i < numImages+1; i++)
  {
    if( i%MOD == 0 ) count++;
    sprintf(textbuffer,"image.%c.%d.gray", (A+i%MOD), (count + OFFSET) );
    printf("%s\n", textbuffer);
    output_file = fopen(textbuffer,"wb");
    fwrite(ptr, 1, IMAGESIZE, output_file);
    fclose(output_file);
    
    ptr += IMAGESIZE + GAPSIZE;
  }
  
  free(filebuffer);
  
  return 0;
}