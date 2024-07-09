// Checkpointing functions

#include <signal.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include "global.h"


static sigset_t set;
static struct sigaction action; 
static int checkpoint;
static int restart_flag;

// Initialization 
void ckpt_register()
{
  sigemptyset(&set);
  sigaddset(&set, SIGURG);
  sigemptyset(&action.sa_mask);
  action.sa_flags = 0;
  action.sa_handler = ckpt_handle_sigurg;
  sigaction(SIGURG, &action, NULL);
}

// Block arrival of signals
void ckpt_block()
{
 sigprocmask(SIG_BLOCK, &set, NULL); 
}

// Allow signals through
void ckpt_unblock()
{
  sigprocmask(SIG_UNBLOCK, &set, NULL); 
}

// Returns 1 if I should checkpoint
int ckpt_check()
{
  return checkpoint;
}

// Clears checkpointing flag
void ckpt_reset()
{
  checkpoint = 0;
  sigaction(SIGURG, &action, NULL);
}

// Custom signal handler
void ckpt_handle_sigurg(int sig)
{
  checkpoint = 1;
  //  signal( sig, SIG_DFL);
}


// Write l values & array(n,n) to restart file
int ckpt_write(int* lvec, double** array, int n, int rank)
{

  FILE* fp; 
  char file[MAXLEN];
  char suffix[5];

  strcpy( file, restart_file);
  //  sprintf( suffix, ".%d", rank);
  sprintf( suffix, "-%d-%d-%d", lvec[0], lvec[1], lvec[2]); 
  strcat( file, suffix);

  if ( (fp = fopen(file, "w")) != NULL ) {
    fwrite(&lvec[0], sizeof(int), 3, fp);
    fwrite(&array[0][0], sizeof(double), n*n, fp);
    printf("Done writing: %d %d %d %d \n", n, lvec[0], lvec[1], lvec[2]); 
    fclose(fp);
    return 1;
  }

  return 0;
} 


// Load l values & array from file
int ckpt_load(int* lvec, double** array, int n, int rank) 
{
  FILE* fp;
  char file[MAXLEN];
  char suffix[5];
  
  strcpy( file, restart_file);
  sprintf( suffix, "-%d-%d-%d", lvec[0], lvec[1], lvec[2]); 
  strcat( file, suffix);
  //  printf("Trying: %s\n", file);
  if ( (fp = fopen(file, "r")) != NULL ) {
    fread(&lvec[0], sizeof(int), 3, fp);
    fread(&array[0][0], sizeof(double), n*n, fp);
    fclose(fp);
    printf("[%d] Ckpt load: %d %d %d \n", rank, lvec[0], lvec[1], lvec[2]); 
    restart_flag = 1;
    return 1;
  }  
  return 0;
}

// Restart inquiry
int ckpt_restart()
{
  return restart_flag;
}
