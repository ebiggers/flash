#include "utilities.h"

/*Merge reads if the best overlap score satisfies the threshold.*/
extern int combine(char id[READ_LENGTH],char read1[READ_LENGTH],char read2[READ_LENGTH],char qual1[READ_LENGTH],char qual2[READ_LENGTH],char read[FINAL_READ_LENGTH],char qual[FINAL_READ_LENGTH],int min_overlap,int max_overlap,float max_mismatch,int phred_offset);

/*Find the best overlap between two reads.*/
extern int pair_align(char read1[READ_LENGTH],char read2[READ_LENGTH],char qual1[READ_LENGTH],char qual2[READ_LENGTH],int min_ovrelap,int max_overlap,int phred_offset,float max_mismatch);

/*Find overlap score for the overlap between two reads at the given position.*/
extern float align_position(int position,char read1[READ_LENGTH],char read2[READ_LENGTH],char qual1[READ_LENGTH],char qual2[READ_LENGTH], float *score,int min_overlap,int max_overlap,int phred_offset);
