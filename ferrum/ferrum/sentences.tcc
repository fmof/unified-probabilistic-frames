#ifndef FERRUM_SENTENCES_TCC_
#define FERRUM_SENTENCES_TCC_

#include <list>
#include <queue>
#include <vector>

namespace ferrum {
  template <typename Predicate>
  DependencyGraph& DependencyGraph::fill(const concrete::DependencyParse& dp, Predicate pred, bool inverse) {
    for(const concrete::Dependency& dep : dp.dependencyList) {
      graph_[dep.gov + 1].push_back(dep.dep + 1);
      edges_[dep.gov + 1][dep.dep+1] = dep.edgeType;
      if(dep.gov < 0) {
	if(! pred(dep.dep) ) continue;
	roots_.emplace(dep.dep + 1);
      }
    }
    if(inverse) {
      // inverse graph
      std::vector< std::list<int> > igraph(num_vertices_, std::list<int>());
      for(const concrete::Dependency& dep : dp.dependencyList) {
	igraph[dep.dep + 1].push_back(dep.gov + 1);
      }
      for(size_t child = 1; child < num_vertices_; ++child) {
	if(! pred(child - 1) ) continue;
	for(int parent : igraph[child]) {
	  if(parent == 0) continue;
	  if(pred(parent - 1)) goto pred_success;
	}
	roots_.emplace(child);
      pred_success:
	continue;
      }
    }
    return *this;
  }

  template <typename ReturnOp>
  void DependencyGraph::bfs_from(int root, ReturnOp& ret_op) const {
    this->bfs_from(root, ret_op, ret_op);
  }
  template <typename ReturnOp, typename NextOp>
  void DependencyGraph::bfs_from(int root, ReturnOp& ret_op, const NextOp& next_op) const {
    std::queue<int> descendants;
    std::set<int> seen;
    descendants.push(root);
    do {
      int curr = descendants.front();
      descendants.pop();
      if(seen.count(curr)) continue;
      seen.emplace(curr);
      ret_op(curr - 1);
      for(int next : graph_[curr]) {
	if(next_op(curr, next, *this)) {
	  descendants.push(next);
	}
      }
    } while(! descendants.empty() );
  }

  template <typename ReturnOp, typename NextOp>
  void DependencyGraph::bfs(ReturnOp& ro, const NextOp& no) const {
    for(int root : roots_) {
      bfs_from(root, ro, no);
    }
  }
}

#endif
