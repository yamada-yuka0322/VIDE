#ifndef __COSMO_QUEUE_HPP
#define __COSMO_QUEUE_HPP

#include <cmath>

namespace CosmoTool {

  template<typename T,typename QType = int>
  class BoundedQueue
  {
  public:
    BoundedQueue(int maxSize, QType defaultMax);
    BoundedQueue(T *pQueue, int maxSize, QType defaultMax);
    ~BoundedQueue();

    void push(T a, QType v);
    T *getQueue() { return m_queue; }
    QType *getPriorities() { return priority; }
    QType getMaxPriority() { return priority[maxQueueSize-1]; }
  private:
    int maxQueueSize;
    T *m_queue;
    QType *priority;
    bool autoFree;
  };
};

#include "bqueue.tcc"

#endif
