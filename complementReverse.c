#include <stdio.h>
#include <string.h>
#include "complementReverse.h"
#include "utilities.h"

/*Get the complement reverse of a string.*/
void compRev(char s[])
{
  int length,i;
  length=strlen(s);

  for(i=0;i<length;i++)
  {
    s[i]=complement(s[i]);
  }
  reverse(s);
}

/*Get bio complement of a character.*/
char complement(char c)
{
  switch (c){
  case 'A':
  case 'a':
    return 'T';
    break;
  case 'C':
  case 'c':
    return 'G';
    break; 
  case 'G':
  case 'g':
    return 'C';
    break;
  case 'T':
  case 't':
    return 'A';
    break;

  default:
  return 'N';
  }
}

/*Get reverse string.*/
void reverse(char s[])
{
  char temp[READ_LENGTH];
  int length,i;

  length=strlen(s);
  for(i=0;i<length;i++){
    temp[i]=s[length-i-1];
  }
  for(i=0;i<length;i++){
    s[i]=temp[i];
  }
  s[length]='\0';
}
