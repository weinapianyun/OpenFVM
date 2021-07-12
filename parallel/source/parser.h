
/* (c) Donald Axel GPL - license */
/* ANSI - C program demonstration, command line calculator */
/* Recursive descent parser */
/* Improve: Make a HELP command. Add more variables.       */

// modified original file calcu.c by OpenFVM team in 12/12/2005

#define MAXL 8196

char *gs;
char *ctp;
char *errorp;
double oldval;
double cx, cy, cz;
double cu, cv, cw, cp, cT, cs;

/* local prototypes(局部原型): */ 
int calcu ();
int evaluate (char *line, double *prev_result);

/* More local prototypes. This could, of course, be a separate file. */
double expression ();
double product ();
double potens ();
double signedfactor ();
double factor ();
double stdfunc ();
