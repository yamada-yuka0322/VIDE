#include <iostream>
#include "bqueue.hpp"

using namespace std;

int main(int argc, char **argv)
{
  CosmoTool::BoundedQueue<int,int> bq(4, 100000.);
  
  for (int i = 10; i >= 0; i--)
    {
      bq.push(i, i);
      
      int *prio = bq.getPriorities();
      for (int j = 0; j < 4; j++)
	cout << prio[j] << " ";

      cout << endl;
    }

  for (int i = 1; i >= -2; i--)
    {
      bq.push(i, i);
      
      int *prio = bq.getPriorities();
      for (int j = 0; j < 4; j++)
	cout << prio[j] << " ";

      cout << endl;
    }

  return 0;
}
