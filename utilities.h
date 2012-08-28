#define READ_LENGTH 170

#define FINAL_READ_LENGTH (2*READ_LENGTH)

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

/*Remove all blanks at the end of the line/string.*/
extern void removeEndBlanks(char s[]); 

/*Convert string to all uppercase letters.*/
extern void uppercase(char s[]); 

/*Set ID of the merged mates.*/
extern void idSet(char id1[READ_LENGTH],char id2[READ_LENGTH],char id[READ_LENGTH]);

/*Print instructions for running the program.*/
extern void printHelp();

/*Print 'help' message.*/
extern void printForHelp();

/*Check that input parameter has a value.*/
extern int checkParity(char parameter, int index, int totalInputs);

/*Calculate the length of read*/
extern int readLength(char read[FINAL_READ_LENGTH]);
