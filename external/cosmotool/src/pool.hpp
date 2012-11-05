#ifndef __COSMO_POOL_HPP
#define __COSMO_POOL_HPP

#include <list>
#include "config.hpp"

namespace CosmoTool
{

  template<typename T>
  struct PoolNode
  {
    T *data;
    uint32_t last_free, size;
    PoolNode<T> *next;
  };

  template<typename T> class MemoryPool;

  template<typename T>
  class MemoryIterator
  {
  private:
    friend  class MemoryPool<T>;

    PoolNode<T> *cur, *previous;
    uint32_t in_node;

    MemoryIterator(PoolNode<T> *h) 
    {
      cur = h; 
      previous = h; 
      in_node = 0;
    }
  public:
    MemoryIterator() { cur = 0; }
    ~MemoryIterator() {}

    const MemoryIterator& operator=(const MemoryIterator& i)
    {
      cur = i.cur;
      previous = i.previous;
      in_node = i.in_node;
    }
    
    bool operator==(const MemoryIterator& i) const
    {
      return (cur == i.cur) && (in_node == i.in_node);
    }

    MemoryIterator& operator++()
    {
      if (cur == 0)
	return *this;

      in_node++;
      if (in_node == cur->size)
	{
	  in_node = 0;
	  previous = cur;
	  cur = cur->next;
	}
      return *this;
    }

    T& operator*()
    {
      return cur->data[in_node];
    }

    T& operator->()
    {
      return cur->data[in_node];
    }

    
    
  };

  // This is bare simple memory pools
  template<typename T>
  class MemoryPool
  {
  private:
    uint32_t m_allocSize;
    PoolNode<T> *head, *current;
    typedef MemoryIterator<T> iterator;
  public:
    MemoryPool(uint32_t allocSize)
      : m_allocSize(allocSize), head(0), current(0) {}

    ~MemoryPool()
    {
      free_all();
    }

    void free_all()
    {
      PoolNode<T> *node = head;

      while (node != 0)
	{
	  PoolNode<T> *next = node->next;
	  	  
	  delete[] node->data;
	  delete node;
	  node = next;
	}
      current = head = 0;
    }

    T *alloc()
    {
      T *ret = alloc_in_node();
      return (ret == 0) ? alloc_new_in_node() : ret;
    }
    
    iterator begin()
    {
      return iterator(head);
    }
    
    iterator end()
    {
      return iterator(0);
    }    

  protected:
    T *alloc_in_node()
    {
      if (current == 0 || current->last_free == current->size)
	return 0;
      return &current->data[current->last_free++];
    }

    T *alloc_new_in_node()
    {
      PoolNode<T> *newNode = new PoolNode<T>;
      if (newNode == 0)
	return 0;

      newNode->last_free = 1;
      newNode->size = m_allocSize;
      newNode->data = new T[m_allocSize];
      if (newNode->data == 0)
	{
	  delete newNode;
	  return 0;
	}
      newNode->next = 0;

      if (current == 0)
	current = head = newNode;
      else
	{
	  current->next = newNode;
	  current = newNode;
	}
      return &newNode->data[0];
    }

  };
  
};

#endif
