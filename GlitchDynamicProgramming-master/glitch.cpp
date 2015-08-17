#include <stdio.h>
#include <iostream>
using namespace std;
int away=0,total=0,ycounter=1,y;
int multiplier=1;
int move_forward(int counter, int steps, int x)
{
    for(int i=1;i<=counter*x;i++)
    {
        if(ycounter%y==0) { multiplier=multiplier*(-1);}
        away+=1*multiplier;
      //  cout<<"\n Away is "<<away;
        ycounter++;
        
    }
    return counter*x;
}
int move_backward(int counter, int steps, int x)
{
    for(int i=1;i<=counter+x;i++)
    {
        if(ycounter%y==0) { multiplier=multiplier*(-1);}
        away-=1*multiplier;
      //  cout<<"\n Away is "<<away;
       ycounter++;
        
    }
    
  
    return counter+x;
}
bool check(int total,int  steps)
{
    if(steps>total) return false;
    else if(steps==total) return false;
    else if(steps<total)
    {
        away-=(steps-total); //not sure
        return true;
    }
    return false;
}
void compute(int x,int steps)
{
    int counter=1;
    while(total<steps)
    {
        total+=move_forward(counter,steps,x);
       // cout<<"\n Away is "<<away;
        if(check(total,steps)) break;
        total+=move_backward(counter,steps,x);
     //   cout<<"\n Away is "<<away;
        if(check(total,steps)) break;
        counter++;
    }

}
int main()
{
    cout<<"Enter x and y for Glitch";
    int x,steps;
    cin>>x>>y;
    cout<<"How many steps has he moved?";
    cin>>steps;
    compute(x,steps);
    cout<<"He is a total of "<<away<<"steps away";
    return 0;
}
