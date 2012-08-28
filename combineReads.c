#include <stdio.h>
#include <string.h>
#include "combineReads.h"
#include "utilities.h"

int combine(char id[READ_LENGTH],char read1[READ_LENGTH],char read2[READ_LENGTH],char qual1[READ_LENGTH],char qual2[READ_LENGTH],char read[FINAL_READ_LENGTH],char qual[FINAL_READ_LENGTH],int min_overlap,int max_overlap,float max_mismatch,int phred_offset){
  int i,j;
  int position; /*starting position of the alignment in the first read, 0-based*/
  int length1,length2; /*lengths of two reads*/
  int overlap; /*length of the overlapping part between two reads*/
  int remain; /*part of the second read not overlapped with the first read*/
  int total; /*length of the final read*/

  position=pair_align(read1,read2,qual1,qual2,min_overlap,max_overlap,phred_offset,max_mismatch);
    if(position<0){
      return 0;
    }
    else{
      length1=strlen(read1);
      length2=strlen(read2);
      overlap=length1-position;
      remain=length2-overlap;
      total=length1+remain;

      /*copy the beginning of the first read*/
      for(i=0;i<position;i++){
	read[i]=read1[i];
	qual[i]=qual1[i];
      }
      /*overlapping part*/
      j=0;
      for(i=position;i<length1;i++){
	if(read1[i]==read2[j]){
	  read[i]=read1[i];
	  qual[i]=max((int)qual1[i],(int)qual2[j]); /*uses ascii code values*/
	}
	else{
	  qual[i]=min((int)qual1[i],(int)qual2[j]);
	  if((int)qual[i]>(2+phred_offset)){
	    qual[i]=2+phred_offset; /*set qual value to phred score of 2*/
	  }
	  if((int)qual1[i]>(int)qual2[j]){ 
	    read[i]=read1[i]; /*take read with higher quality value*/
	  }
	  else{
	    if((int)qual1[i]<(int)qual2[j]){
	      read[i]=read2[j];
	    }
	    else{ /*same quality values--> check if any base is N*/
	      if(read2[j]=='N'){
		read[i]=read1[i];
	      }
	      else{
		read[i]=read2[j];
	      }
	    }
	  }
	}
	j++;
      }
      /*copy the remainder of the second read*/
      for(i=length1;i<total;i++){
	read[i]=read2[j];
	qual[i]=qual2[j];
	j++;
      }
      read[total]='\0';
      qual[total]='\0';
    }
    return 1;
}

/*Return the position of the first overlapping character in the first read (0-based) of the best overlap between two reads, or -1 if the best overlap does not satisfy the required threshold.*/
int pair_align(char read1[],char read2[],char qual1[],char qual2[],int min_overlap,int max_overlap,int phred_offset,float max_mismatch){
  int i,read1_len,best_position;
  float qual_score,best_qual_score,score,best_mismatch;

  read1_len=strlen(read1);
  best_qual_score=align_position(0,read1,read2,qual1,qual2,&score,min_overlap,max_overlap,phred_offset);
  if(qual_score>-1){
    best_mismatch=score;
    best_position=0;
  }
  else{
    best_mismatch=10000; /*default: bad value*/
    best_position=-1;
  }

  for(i=1;i<(read1_len-min_overlap+1);i++){ /* requires at least min_overlap bases overlap*/
    qual_score=align_position(i,read1,read2,qual1,qual2,&score,min_overlap,max_overlap,phred_offset);
    if(score<best_mismatch){ 
      best_qual_score=qual_score;
      best_mismatch=score;
      best_position=i; 
    }
    else{
      if((score==best_mismatch) && (qual_score<best_qual_score)){
        best_qual_score=qual_score;
        best_position=i;
      }
    }
  }

  if(best_mismatch>max_mismatch){
    return -1;
  }
  else{
    return best_position;
  }
}

/*Align the read2 to the read1 starting at the given position (0-based) in the read1. Return score=mismatches/overlap_length and qual_score=average of quality scores of all mismatches in the overlapping region. Max overlap length considered in the scoring is max_ovrelap*/
float align_position(int position,char read1[],char read2[],char qual1[],char qual2[],float *score,int min_overlap,int max_overlap,int phred_offset){
  int j,read1_len,read2_len,overlap,not_overlap;
  float qual_score=0.0; /*qual_score is for mismatches only, so lower qual_score is preferable since lower qual_score means less confidence in base calls at mismatches*/
  *score=0;
  read1_len=strlen(read1);
  overlap=read1_len-position;
  read2_len=strlen(read2);
  if(read2_len<overlap){ /*read2 is contained in read1->read1 can't be extended->return -1, so extension is not considered*/
    *score=10001;
    return -1.0;
  }

  not_overlap=0;
  for(j=0;j<overlap;j++){
    if((read1[j+position]=='N') || (read1[j+position]=='n') || (read2[j]=='N') || (read2[j]=='n')){
      not_overlap++;
    }
    else{
      if(read1[j+position]!=read2[j]){
	*score=*score+1;
	qual_score=qual_score+min(qual1[j+position],qual2[j])-phred_offset;
      }
    }
  }
  overlap-=not_overlap; /*reduce the length of overlap by the number of matching N's*/
  if(overlap>max_overlap){
    qual_score=qual_score/max_overlap;
    *score=*score/max_overlap;
  }
  else{
    if(overlap>=min_overlap){
      qual_score=qual_score/overlap;
      *score=*score/overlap;
    }
    else{
      *score=1000000; /*higher value than in default when position=0, so that it does not change best_position*/
      return -1;
    }
  }
  return qual_score;
}
