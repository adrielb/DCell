/*
 * Spool Header:  2*2752 Bytes = 5504 Bytes
 * Image Gap:     2*128 Bytes = 256 Bytes
 * Image Size:    2*512*512 = 524288 Bytes
 * 
 */

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#define HEADERSIZE 5504
#define IMAGESIZE 2*256*256
#define GAPSIZE 256
#define TEXTSIZE 1024

int main(int argc, char **argv) 
{
	char *rootdir = getcwd( NULL, 0 );
	FILE *input_file, *output_file;
	char *spoolfilename = argv[1];
	
	unsigned short image[IMAGESIZE];
	char gap[GAPSIZE];
	char textbuffer[TEXTSIZE];
	char *filebuffer;
	
	printf("Processing: %s\n", spoolfilename);
		
	if( argc != 2 )
	{
		fprintf(stderr, "Requires spool file name\n");
		free(rootdir);
		return 1;
	}
		
	input_file = fopen(spoolfilename,"rb");
  
	if (input_file == NULL) {
		fprintf(stderr, "Cannot open: %s\n", spoolfilename);
		free(rootdir);
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
  	fprintf(stderr, "Size read (%lu) does not equal to file size (%lu)",
  		s, mystat.st_size);
  	return 3;
  }
  
  /* Create directory based on spool file name*/
  sprintf(textbuffer, "DIR.%s", spoolfilename);
	mkdir(textbuffer, S_IRWXU|S_IRWXG);
	chdir(textbuffer);
	rootdir = getcwd( NULL, 0);
  
  output_file = fopen("header","wb");
	fwrite(filebuffer, 1, HEADERSIZE, output_file);
	fclose(output_file);
	
	const int OFFSET = 10000;
	int numImages = (mystat.st_size - HEADERSIZE) / (GAPSIZE + IMAGESIZE);
	char *ptr = filebuffer + HEADERSIZE;
	for( int count = 1; count < numImages / 1; count++)
	{  	
//	printf("%lu\n",mystat.st_size - (ptr - filebuffer));
		sprintf(textbuffer,"image.a.%d.gray", (count+OFFSET) );
		output_file = fopen(textbuffer,"wb");
		fwrite(ptr, 1, IMAGESIZE, output_file);
		fclose(output_file);
		printf("%s\n",textbuffer);
		
		ptr += IMAGESIZE + GAPSIZE;
/*		
		sprintf(textbuffer,"image.b.%d.gray", (count+OFFSET) );
		output_file = fopen(textbuffer,"wb");
		fwrite(ptr, 1, IMAGESIZE, output_file);
		fclose(output_file);
		printf("%s\n",textbuffer);
		
		ptr += IMAGESIZE + GAPSIZE;
*/
	}
	
	free(filebuffer);
	chdir("..");
	
	return 0;
	
	for( int count = 1; count < numImages/2; count++)
	{
		fread(image, 2, IMAGESIZE, input_file);
		
		sprintf(textbuffer,"image.a.%d.gray",count );
		output_file = fopen(textbuffer,"wb");
		fwrite(image, 1, IMAGESIZE, output_file);
		fclose(output_file);
		printf("%s\n",textbuffer);
		fread(gap, 1, GAPSIZE, input_file);
		
		
		/*fread(image, 2, IMAGESIZE, input_file);
		
		sprintf(buffer,"%simage.b.%d.gray",rootdir,count );
		output_file = fopen(buffer,"wb");
		fwrite(image, 2, IMAGESIZE, output_file);
		fclose(output_file);
		printf("%s\n",buffer);
		fread(gap, 1, GAPSIZE, input_file);
		*/
		
	}
	return 0;	
}
