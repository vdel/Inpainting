#ifndef OUTPUTH
#define OUTPUTH

#include <math.h>

class ProgressBar
{
  private:
  int min, max;
  unsigned int last_pct;
  
  void write(int pct)
  {
    last_pct = pct;
    printf("%d%%",pct);
    fflush(stdout);       
  }
  
  void backtrack()
  {    
    int l = ((last_pct==0)?0:int(log10(last_pct)))+1;
    printf("\033[%dD", l+1);
  }

  public:  
  ProgressBar(int _min, int _max)
  {
    min = _min;
    max = _max;
    write(0);
  }  

  void update(int value)
  {
    unsigned int pct = (value-min)*100/(max-min);
    if(pct != last_pct)
    {
      backtrack();
      write(pct);
    }
  }
  
  void done()
  {
    backtrack();
    printf("done.\n");
    fflush(stdout);
  }
};

#endif
