#ifndef SIMPLE_BVH_H
#define SIMPLE_BVH_H
#include "clipper3/clipper.h"
#include <vector>
#include <numeric>

struct BVHNode {
    Clipper2Lib::Rect64 box;
    int left{-1};
    int right{-1};
    int index{-1};
};

struct BVH {
    std::vector<BVHNode> nodes;
};

inline Clipper2Lib::Rect64 combine(const Clipper2Lib::Rect64& a,const Clipper2Lib::Rect64& b){
    return {std::min(a.left,b.left), std::min(a.bottom,b.bottom),
            std::max(a.right,b.right), std::max(a.top,b.top)};
}

inline BVH buildBVH(const Clipper2Lib::Paths64& paths){
    BVH bvh;
    size_t n=paths.size();
    if(n==0) return bvh;
    bvh.nodes.reserve(2*n);
    std::vector<int> indices(n);
    std::iota(indices.begin(),indices.end(),0);
    std::function<int(int,int)> build=[&](int l,int r){
        int nodeIdx=bvh.nodes.size();
        bvh.nodes.push_back({});
        if(l==r){
            bvh.nodes[nodeIdx].index=indices[l];
            bvh.nodes[nodeIdx].box=Clipper2Lib::Rect64{0,0,0,0};
            if(!paths[indices[l]].empty()) bvh.nodes[nodeIdx].box=Clipper2Lib::GetBounds(paths[indices[l]]);
            return nodeIdx;
        }
        int mid=(l+r)/2;
        int left=build(l,mid);
        int right=build(mid+1,r);
        bvh.nodes[nodeIdx].left=left;
        bvh.nodes[nodeIdx].right=right;
        auto bl=bvh.nodes[left].box;
        auto br=bvh.nodes[right].box;
        bvh.nodes[nodeIdx].box=combine(bl,br);
        return nodeIdx;
    };
    build(0,n-1);
    return bvh;
}

inline bool overlapBVH(const BVH& A,const BVH& B,
                       const Clipper2Lib::Paths64& pA,
                       const Clipper2Lib::Paths64& pB,int aIdx,int bIdx){
    const BVHNode& na=A.nodes[aIdx];
    const BVHNode& nb=B.nodes[bIdx];
    if(na.box.right < nb.box.left || na.box.left > nb.box.right ||
       na.box.top < nb.box.bottom || na.box.bottom > nb.box.top)
        return false;
    if(na.index!=-1 && nb.index!=-1){
        Clipper2Lib::Paths64 isect=Clipper2Lib::Intersect({pA[na.index]},{pB[nb.index]},Clipper2Lib::FillRule::NonZero);
        return std::abs(Clipper2Lib::Area(isect))>0;
    }
    if(na.index==-1){
        if(overlapBVH(A,B,pA,pB,na.left,bIdx)) return true;
        if(overlapBVH(A,B,pA,pB,na.right,bIdx)) return true;
    } else {
        if(overlapBVH(A,B,pA,pB,aIdx,nb.left)) return true;
        if(overlapBVH(A,B,pA,pB,aIdx,nb.right)) return true;
    }
    return false;
}

inline bool overlapBVH(const BVH& A,const BVH& B,
                       const Clipper2Lib::Paths64& pA,
                       const Clipper2Lib::Paths64& pB){
    if(A.nodes.empty()||B.nodes.empty()) return false;
    return overlapBVH(A,B,pA,pB,0,0);
}

#endif
