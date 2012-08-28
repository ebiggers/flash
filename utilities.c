#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "utilities.h"

/*Remove all blank characters at the end of the given string*/
void removeEndBlanks(char s[READ_LENGTH]){
  int last,i;
  last=strlen(s)-1;
  char c;
  for(i=last;i>=0;i--){
    if((s[i]!='\n') && (s[i]!='\t') && (s[i]!=' ')){
      last=i+1; /*position for end-of-string character*/
      i=-2;
    } 
  }
  if(i==-1){  /*The entire line consists only of blanks*/
    last=0;
  }
  strncpy(s,s,last);
  s[last]='\0';
}

/*Convert string to all uppercase letters.*/
void uppercase(char c[READ_LENGTH]){
  int i,length;
  length=strlen(c);
  for(i=0;i<length;i++){
    c[i]=toupper(c[i]);
  }
}

/*Set up id of the merged mates.*/
void idSet(char id1[READ_LENGTH],char id2[READ_LENGTH],char id[READ_LENGTH]){
  int i,j,len,compareStr,diff;
  char temp[READ_LENGTH];
  compareStr=strcmp(id1,id2);
  if(compareStr==0){
      strcpy(id,id1);
  }
  else{
    i=0;
    while((id1[i]!='/') && (id1[i]!='\n') && (id1[i]!='\0')){
      i++;
    }
    strncpy(id,id1,i);
    id[i]='\0';
    len=strlen(id1);
    j=i+2;
    if(len>j){
      if(id1[j]=='#'){  /*read id has a barcode*/      
	i=j;
 	while((id1[j]!='\n') && (id1[j]!='\0')){
	  j++;
	}
	diff=j-i;
	strncpy(temp,id1+i,diff);
	temp[diff]='\0';
	strcat(id,temp);
	i=j-2;
      }
    }
    id[i]='\0';
  }
}

/*Print instructions for running the program.*/
void printHelp(){
  printf("\n");
  printf("USAGE: flash <mates1.fastq> <mates2.fastq> [options]\n\n");
  printf("DESCRIPTION:\n");
  printf("FLASH (Fast Length Adjustment of SHort reads) is a very accurate and fast tool \n to merge paired-end reads that were generated from DNA fragments whose lengths \n are shorter than twice the length of reads. Merged read pairs result in unpaired \n longer reads. Longer reads are generally more desired in genome assembly and \n genome analysis processes.\n\n");
  printf("MANDATORY INPUT:\n");
  printf("Two fastq files of paired-end reads. Corresponding pairs must be in the same order in both files.\n\n");
  printf("OPTIONS:\n");
  printf("-m: minOverlap is the minimum required overlap length between two reads to provide \n a confident overlap. Default: 10bp.\n");
  printf("-M: maxOverlap is the maximum overlap length expected in approximately 90%% of read \n pairs. It is by default set to 70bp, which works well for 100bp reads generated from \n 180bp library (normal distribution of fragment lengths is assumed). Overlaps longer \n than maxOverlap are still considered as good overlaps, but the mismatch ratio \n (explained below) is calculated over the maxOverlap rather than the true overlap \n length. If you enter a value for maxOverlap, read length, fragment length, and fragment \n standard deviation that you enter will be ignored for calculation of maxOverlap parameter. \n Default: 70bp. \n");
  printf("-x: mismatchRatio is the maximum allowed ratio of the number of mismatches and the \n overlap length. An overlap with mismatch ratio higher than the set value is considered \n incorrect overlap and mates will not be merged. Any occurence of an 'N' in any read is \n ignored and not counted towards the mismatches or overlap length. Our experimental \n results suggest that higher values of mismatchRatio yield larger number of correctly \n merged read pairs but at the expense of higher number of incorrectly merged read \n pairs. Default: 0.25. \n");
  printf("-p: phredOffset is the smallest ASCII value of the characters used to represent \n quality values of bases in fastq files. It should be set to either 33, which corresponds \n to the later Illumina platforms and Sanger platforms, or 64, which corresponds to \n the earlier Illumina platforms. Default: 33.\n");
  printf("-o prefix of output files. Default: out.\n");
  printf("-d path to directory for output files. Default: current working directory.\n");
  printf("-r average reads length. Default: 100.\n");
  printf("-f average fragment length. Default: 180.\n");
  printf("-s standard deviation of fragment lengths. If you do not know standard deviation of the fragment library,\n you can probably assume that the standard deviation is 10%% of the average fragment length. \n Default: 20. \n");
  printf("-h: display help and exit.\n\n");
}

/*Print help message*/
void printForHelp(){
  printf("Try 'flash -h' for more information.\n");
}

/*Check that each input parameter has a value*/
int checkParity(char parameter, int index, int totalInputs){
  if(index>=totalInputs){
    printf("Value of the parameter %c is missing.\n",parameter);
    printForHelp();
    return(1);
  }
  else{
    return(0);
  }
}

/*Calculate read length*/
int readLength(char read[FINAL_READ_LENGTH]){
  int len;
  len=strlen(read);
  len=len-1;
  return len;
}
