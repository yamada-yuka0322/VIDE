#include <list>
#include <queue>

namespace CosmoTool
{

  namespace simple_bsp
  {
    
    template<typename T, typename PType, int N>
    BSP<T,PType,N>::BSP()
     throw()
    {
      root = 0;
    }

    template<typename T, typename PType, int N>
    BSP<T,PType,N>::~BSP()
    {
      while (!allocated.empty())
	{
	  node_t *r = allocated.front();
	      
	  allocated.pop();
	  delete r;
	}
    }

    template<typename T, typename PType, int N>
    void BSP<T,PType,N>::insert(Facet<space_t>& f)
    {
      std::list<node_t **> subtrees;
      std::list<node_t **> leaf_insert;

      if (root != 0)
	subtrees.push_back(&root);
      else
	leaf_insert.push_back(&root);

      // Find the point of insertion. Do not bother to split triangle for this
      // implementation.
      while (!subtrees.empty())
	{
	  std::list<node_t **> new_subtrees;
	  typename std::list<node_t **>::iterator iter = subtrees.begin();

	  while (iter != subtrees.end())
	    {
	      typename space_t::point_t dp;
	      bool cond_plus = false, cond_minus = false;
	      node_t *current = *(*iter);
	      
	      for (int j = 0; j < N; j++)
		{
		  dp = dot_product<space_t>(f.p[j],current->plane.n) + current->plane.d;
		  cond_plus = cond_plus || (dp > 0);
		  cond_minus = cond_minus || (dp <= 0);
		}
	      
	      bool joint = cond_plus && cond_minus;
	      bool not_joint = (!cond_plus) && (!cond_minus);
	      
	      if (joint || not_joint)
		{
		  // Crawl and add another subtree
		  *iter = &(current->minus);
		  if (current->plus != 0)
		    new_subtrees.push_back(&current->plus);
		  else
		    leaf_insert.push_back(&current->plus);
		}
	      else
		{
		  if (cond_plus)
		    *iter = &current->plus;
		  else
		    *iter = &current->minus;
		}
	      if (*(*iter) == 0)
		{
		  leaf_insert.push_back(*iter);
		  iter = subtrees.erase(iter);
		}
	      else
		++iter;
	    }
	  if (!new_subtrees.empty())
	    subtrees.splice(subtrees.end(), new_subtrees);
	}
      
      node_t * current = new node_t;
      f.normal(current->plane.n);
      current->plane.d = -dot_product<space_t>((f.p[0]),current->plane.n);
      
      for (typename std::list<node_t **>::iterator i = leaf_insert.begin();
	   i != leaf_insert.end();
	   ++i)
	*(*i) = current;

      allocated.push(current);
    }

    template<typename T, typename PType, int N>
    bool BSP<T,PType,N>::inside(const typename space_t::coord_t& p) const
    {
      node_t *current = root;

      do
	{
	}
      while();
      current
    }

  };

};
