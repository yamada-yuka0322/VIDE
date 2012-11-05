namespace CosmoTool
{
  
  template<typename T, typename QType>
  BoundedQueue<T,QType>::BoundedQueue(int maxSize, QType defaultVal)
  {
    maxQueueSize = maxSize;
    m_queue = new T[maxSize];
    priority = new QType[maxSize];
    autoFree = true;

    for (int i = 0; i < maxSize; i++)
      priority[i] = defaultVal;
  }
 
  template<typename T, typename QType>
  BoundedQueue<T,QType>::BoundedQueue(T *pQueue, int maxSize, QType defaultVal)
  {
    maxQueueSize = maxSize;
    m_queue = pQueue;
    priority = new QType[maxSize];
    autoFree = false;
    
    for (int i = 0; i < maxSize; i++)
      priority[i] = defaultVal;

    
  }
  
  template<typename T, typename QType>
  BoundedQueue<T,QType>::~BoundedQueue()
  {
    if (autoFree)
      delete[] m_queue;
    delete[] priority;
  }


  template<typename T, typename QType>
  void BoundedQueue<T,QType>::push(T a, QType v)
  {
    if (v > priority[maxQueueSize-1])
      return;

    int i;
    for (i = maxQueueSize-2; i >= 0; i--)
      {
	if (v > priority[i])
	  {
	    priority[i+1] = v;
	    m_queue[i+1] = a;
	    return;
	  }
	else
	  {
	    priority[i+1] = priority[i];
	    m_queue[i+1] = m_queue[i];
	  }
      }
    if (i < 0)
      {
	priority[0] = v;
	m_queue[0] = a;
      }
  }
   
};
