#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "complementReverse.h"
#include "combineReads.h"
#include "utilities.h"
#include <zlib.h>

/*Usage: flash <mates1.fastq> <mates2.fastq> [-m minOverlap] [-M maxOverlap] [-x mismatchRatio] [-p phredOffset] [-o prefixOfOutputFiles] [-d pathToDirectory] [-r readLength] [-f fragmentlength] [-s fragmentStandardDeviation]*/  
/*Try 'flash -h' for more information.*/

int main(int argc, char **argv)
{ 
  int result,i,compareStr,checkPar;
  int start,j,k,dir,stop,length,position,max_length,max_hist;
  int min_overlap,max_overlap,max_overlap_flag,phred_offset,read_length,frag_length;
  int len[FINAL_READ_LENGTH];
  float max_mismatch, frag_std;
  char prefix[30]="out";
  char read[FINAL_READ_LENGTH];
  char outFile[60],notCombined1[60],notCombined2[60],hist[60],histogram[60];
  char directory[200]=".";
  char temp[20][200];
  char tempPath[200]="";
  char blank[1]="";
  char id[READ_LENGTH];
  char file1[4][READ_LENGTH],file2[4][READ_LENGTH],finalFile[4][FINAL_READ_LENGTH];
  FILE *outFile_p,*notCombined1_p,*notCombined2_p,*hist_p,*histogram_p;
  gzFile mates1_g,mates2_g;

  struct stat st;

  phred_offset=33;
  max_overlap=70; /*max overlap value used to calculate mismatch score*/
  max_overlap_flag=0;
  max_mismatch=0.25; /*max ratio=(number of mismatches)/(overlap length)*/
  min_overlap=10;
  read_length=100;
  frag_length=180;
  frag_std=20;

  if(argc==2){
    if((argv[1][0]=='-') && (argv[1][1]=='h')){
      printHelp();
      return(1);
    }
    else{
      printf("At least one input file is not specified.\n");
      printForHelp();
      return(1);
    }
  }
  else{
    if(argc<3){
      printf("At least one input file is not specified.\n");
      printForHelp();
      return(1);
    }
  }

  if(argv[1][0]=='-'){
    printf("One of two input files is not specified or its name starts with '-'. The names of these files cannot start with '-'.\n");
    printForHelp();
    return(1);
  }
  mates1_g=gzopen(argv[1],"r");
  if(mates1_g == NULL){
    printf("The first specified file cannot be opened, or is empty.\n");
    printForHelp();
    return(1);
  }
  if(argv[2][0]=='-'){
    printf("One of two input files is not specified or its name starts with '-'. The names of these files cannot start with '-'.\n");
    printForHelp();
    return(1);
  }
  mates2_g=gzopen(argv[2],"r");
  if(mates2_g == NULL){
    printf("The second specified file cannot be opened, or is empty.\n");
    printForHelp();
    return(1);
  }

  for(i=3;i<argc;i++){
    if(argv[i][0]=='-'){
      switch(argv[i][1]){
        case 'm':
          i++;
	  checkPar=checkParity('m',i,argc);
	  if(checkPar==1){
	    return(1);
	  }
	  else{
	    min_overlap=atoi(argv[i]);
	    if(min_overlap<1){
	      printf("Minimum overlap should be positive.\n");
	      printForHelp();
	      return(1);
	    }
	    break;
	  }
        case 'M':
          i++;
	  checkPar=checkParity('M',i,argc);
	  if(checkPar==1){
	    return(1);
	  }
	  else{
	    max_overlap=atoi(argv[i]);
	    max_overlap_flag=1;
	    if(max_overlap<1){
	      printf("Maximum overlap should be positive.\n");
	      printForHelp();
	      return(1);
	    }
	    break;
	  }
        case 'x':
          i++;
	  checkPar=checkParity('x',i,argc);
	  if(checkPar==1){
	    return(1);
	  }
	  else{
	    max_mismatch=atof(argv[i]);
	    if((max_mismatch<0) || (max_mismatch>1)){
	      printf("Max mismatch rate has to be in the interval [0,1].\n");
	      printForHelp();
	      return(1);
	    }
	    break;
	  }
        case 'o':
	  i++;
	  checkPar-checkParity('o',i,argc);
	  if(checkPar==1){
	    return(1);
	  }
	  else{
	    strcpy(prefix,argv[i]);
	    break;
	  }
        case 'd':
	  i++;
	  checkPar=checkParity('d',i,argc);
	  if(checkPar==1){
	    return(1);
	  }
	  else{
	    strcpy(directory,argv[i]);
	    break;
	  }
        case 'p':
          i++;
	  checkPar=checkParity('p',i,argc);
	  if(checkPar==1){
	    return(1);
	  }
	  else{
	    phred_offset=atoi(argv[i]);
	    if((phred_offset!=33) && (phred_offset!=64)){
	      printf("WARNING: Phred offset is usually either 64 (for earlier Illumina data) or 33 (for Sanger and later Illumina data).\n"); 
	      printForHelp();
	    }
	    break;
	  }
        case 'r':
	  i++;
	  checkPar=checkParity('r',i,argc);
	  if(checkPar==1){
	    return(1);
	  }
	  else{
	    read_length=atoi(argv[i]);
	    if((read_length<=0) || (read_length>READ_LENGTH)){
	      printf("Read length has to be in the range (0,170]. If the reads are longer than 170, please change the constant value for READ_LENGTH in utilities.h\n");
	      printForHelp();
	      return(1);
	    }
	    break;
	  }
        case 'f':
	  i++;
	  checkPar=checkParity('f',i,argc);
	  if(checkPar==1){
	    return(1);
	  }
	  else{
	    frag_length=atoi(argv[i]);
	    if(frag_length<=0){
	      printf("Fragment length must be positive.\n");
	      printForHelp();
	      return(1);
	    }
	    break;
	  }
        case 's':
	  i++;
	  checkPar=checkParity('s',i,argc);
	  if(checkPar==1){
	    return(1);
	  }
	  else{
	    frag_std=atof(argv[i]);
	    if(frag_std<=0){
	      printf("Standard deviation of fragments must be positive.\n");
	      printForHelp();
	      return(1);
	    }
	    break;
	  }
        case 'h':
	  printHelp();
	  return(1);
	  break;
        default:
          printf("Unknown option. %s\n.",argv[i]);
          printf("Usage: flash <mates1.fastq> <mates2.fastq> [-m minOverlap] [-M maxOverlap] [-x mismatchRatio] [-p phredOffset] [-o prefixOfOutputFiles] [-d pathToOutputDirectory] [-r readLength] [-f fragmentLength] [-s fragmentStandardDeviation]\n");
	  printForHelp();
          return(1);
      }
    }
    else{  /*else of if(argv[i][0]=='-'*/
      printf("Unknown input parameter %s.\n", argv[i]);
      printForHelp();
    }
  }

  if(max_overlap_flag==0){
    if((read_length!=100) || (frag_length!=180) || (frag_std!=20)){
      max_overlap=(int)(2*read_length-frag_length+2.5*frag_std);
    }
  }

  if(directory!="."){
    j=0;
    strcpy(temp[0],directory);
    start=0;
    dir=strlen(directory);
    for(i=0;i<dir;i++){
      if((directory[i]=='/') | (directory[i]=='\\')){
	strncpy(temp[j],directory+start,i-start);
	temp[j][i]='\0';
	j++;
	start=i+1;
      }
    }
    strncpy(temp[j],directory+start,i-start);
    temp[j][i]='\0';
    j++;
    stop=0;

    /*check if the path is given as absolute path*/
    if((directory[0]=='/') || (directory[0]=='\\')) { /*absolute path is provided*/
      /*check which part of absolute path exists*/
      i=1;
      strcpy(tempPath,blank);
      strcat(tempPath,"/");
      strcat(tempPath,temp[1]);
      while(stat(tempPath,&st)==0){
	i++;
	strcat(tempPath,"/");
	strcat(tempPath,temp[i]);
      }
      stop=i;
      i=j;
    }

      /*change directory to the latest existing stage of absolute path*/
    if(stop==1){
      strcpy(tempPath,"/");
      printf("WARNING: You are trying to create a directory at the root. If you receive a 'Segmentatoin fault' message, you most probably do not have right to create a directory at root. \n");
    }
    else{
      strcpy(tempPath,blank);
      for(i=1;i<stop;i++){
	strcat(tempPath,"/");
	strcat(tempPath,temp[i]);
      }
    }
    chdir(tempPath);

    /*create directories for output*/
    for(i=stop;i<j;i++){
      mkdir(temp[i], S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      chdir(temp[i]);
    }
  }

  strcpy(notCombined1,prefix);
  strcat(notCombined1,".notCombined_1.fastq");
  strcpy(notCombined2,prefix);
  strcat(notCombined2,".notCombined_2.fastq");
  strcpy(outFile,prefix);
  strcat(outFile,".extendedFrags.fastq");
  notCombined1_p=fopen(notCombined1,"w");
  notCombined2_p=fopen(notCombined2,"w");
  outFile_p=fopen(outFile,"w");

  gzgets(mates1_g,file1[0],READ_LENGTH);
  removeEndBlanks(file1[0]);
  while(!gzeof(mates1_g)){
    for(i=1;i<4;i++){
      gzgets(mates1_g,file1[i],READ_LENGTH);
      removeEndBlanks(file1[i]);
    }
    for(i=0;i<4;i++){
      gzgets(mates2_g,file2[i],READ_LENGTH);
      removeEndBlanks(file2[i]);
    }
    uppercase(file1[1]);
    compRev(file2[1]);
    reverse(file2[3]);
    idSet(file1[0],file2[0],id);
    result=combine(id,file1[1],file2[1],file1[3],file2[3],finalFile[1],finalFile[3],min_overlap,max_overlap,max_mismatch,phred_offset);
    if(result){
      fprintf(outFile_p,"%s\n%s\n+\n%s\n",id,finalFile[1],finalFile[3]);
    }
    else{
      fprintf(notCombined1_p,"%s\n%s\n+\n%s\n",file1[0],file1[1],file1[3]);
      compRev(file2[1]);
      reverse(file2[3]);
      fprintf(notCombined2_p,"%s\n%s\n+\n%s\n",file2[0],file2[1],file2[3]);
    }
    gzgets(mates1_g,file1[0],READ_LENGTH); /*get id of the next read. Do it here because if it is EOF, it will not go into while-loop*/
      removeEndBlanks(file1[0]);
  }

  gzclose(mates1_g);
  gzclose(mates2_g);
  fclose(outFile_p);
  fclose(notCombined1_p);
  fclose(notCombined2_p);

  /*Make histogram of lengths of merged reads*/
  for(i=read_length;i<(FINAL_READ_LENGTH+1);i++){
    len[i]=0;
  }
  max_length=0;
  outFile_p=fopen(outFile,"r");
  fgets(read,FINAL_READ_LENGTH,outFile_p);
  while(!feof(outFile_p)){
      fgets(read,FINAL_READ_LENGTH,outFile_p);
      length=readLength(read);
    if(length>max_length){
      max_length=length;
    }
    len[length]++;
    for(i=0;i<3;i++){
      fgets(read,FINAL_READ_LENGTH,outFile_p);
    }
  }
  fclose(outFile_p);

  max_hist=0;
  strcpy(hist,prefix);
  strcat(hist,".hist");
  hist_p=fopen(hist,"w");
  for(i=read_length;i<=max_length;i++){
    fprintf(hist_p,"%d\t%d\n",i,len[i]);
    if(len[i]>max_hist){
      max_hist=len[i];
    }
  }
  fclose(hist_p);

  j=max_hist/100;
  if(j==0){
    j=1;
  }
  strcpy(histogram,prefix);
  strcat(histogram,".histogram");
  histogram_p=fopen(histogram,"w");
  hist_p=fopen(hist,"r");
  fscanf(hist_p,"%d %d",&position,&length);
  while(!feof(hist_p)){
    if(length==0){
      k=0;
    }
    else{
      k=length/j;
      if(k==0){
	k=1;
      }
    }
    fprintf(histogram_p,"%d\t",position);
    for(i=0;i<k;i++){
      fprintf(histogram_p,"*");
    }
    fprintf(histogram_p,"\n");
    fscanf(hist_p,"%d %d",&position,&length);
  }
  fclose(hist_p);
  fclose(histogram_p);

  return(0);
}
