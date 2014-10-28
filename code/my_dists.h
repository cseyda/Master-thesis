#ifndef MYDISTS_H
#define MYDISTS_H

#include <cmath>
#include <cstdlib>
#include <string.h>
#include <stdint.h>

#include <iostream>

#include "flann/flann.hpp"

using namespace std;
using namespace flann;

template<class T>
struct Cosine
{
    //typedef bool is_vector_space_distance;

    typedef T ElementType;
    typedef typename Accumulator<T>::Type ResultType;

    template <typename Iterator1, typename Iterator2>
    ResultType operator()(Iterator1 a, Iterator2 b, size_t size, ResultType /*worst_dist*/ = -1) const
    {
        ResultType ab = ResultType();
        ResultType a_sq = ResultType();
        ResultType b_sq = ResultType();
        
        
        ResultType a_;
        ResultType b_;
        for(size_t i = 0; i < size; ++i ) {
            a_ = *a++;
            b_ = *b++;
            ab += a_ * b_;
            a_sq += a_ * a_;
            b_sq += b_ * b_;
        }
        /*
        Iterator1 last = a + size;
        Iterator1 lastgroup = last - 3;
*/
        /* Process 4 items with each loop for efficiency. 
        while (a < lastgroup) {
            //diff0 = (ResultType)std::abs(a[0] - b[0]);
            ab += (ResultType) a[0] * b[0] 
                + a[1] * b[1] 
                + a[2] * b[2] 
                + a[3] * b[3];
            
            a_sq += (ResultType) a[0] * a[0] 
                + a[1] * a[1] 
                + a[2] * a[2] 
                + a[3] * a[3];
            
            b_sq += (ResultType) b[0] * b[0] 
                + b[1] * b[1] 
                + b[2] * b[2] 
                + b[3] * b[3];
            
            a += 4;
            b += 4;
        }
        
        ResultType a_;
        ResultType b_;*/
        /* Process last 0-3 pixels.  Not needed for standard vector lengths. 
        while (a < last) {
            a_ = *a++;
            b_ = *b++;
            ab += (ResultType) a_ * b_;
            a_sq += (ResultType) a_ * a_;
            b_sq += (ResultType) b_ * b_;
        }*/
        
        ResultType res = (ResultType) (1.0 - (ab / (sqrt(a_sq) * sqrt(b_sq))));
        //std::cout << res << " " << ab << " " << a_sq << " " << b_sq << std::endl;
        return res;
    }
    /*
    template <typename U, typename V>
    inline ResultType accum_dist(const U& a, const V& b, int) const
    {
        return (a-b)*(a-b);
    }
    */
};

template<class T>
struct Jaccard
{
    typedef bool is_vector_space_distance;

    typedef T ElementType;
    typedef typename Accumulator<T>::Type ResultType;

    template <typename Iterator1, typename Iterator2>
    //ResultType operator()(Iterator1 a, Iterator2 b, size_t size, ResultType /*worst_dist*/ = -1) const
    /*{
        
        unsigned int f_01 = 0;
        unsigned int f_11 = 0; 
                
        for(size_t i = 0; i < size; ++i ) {
            if (*a++ > 0) {
                if (*b++ > 0){
                    f_11++;
                } else {
                    f_01++;
                }
            } else {
                if (*b++ > 0){
                    f_01++;
                }
            }
        }
        
        ResultType res = f_01+f_11;
        if ( res != 0) {
            res = (ResultType) (f_01) / res;
        } else {
            res = 1;
        }
        
        return res;
    };*/
    
    ResultType operator()(Iterator1 a, Iterator2 b, size_t size, ResultType /*worst_dist*/ = -1) const
    /**/{
        
        unsigned int f_01 = 0;
        unsigned int f_11 = 0; 
        
        //for(size_t i = 0; i < size; ++i ) {
        while((*a != 0) && (*b !=0)){
            if (*a == 0){
                while (*b != 0){
                    b++;
                    f_01++;
                }
                break;
            } else if ((*b) == 0){
                while (*a != 0){
                    a++;
                    f_01++;
                }
                break;
            }
            
            else if (*a > *b) {
                b++;
                f_01++;
            } else if (*a < *b) {
                a++;
                f_01++;
            } else {//found pair
                a++;
                b++;
                f_11++;
            }
        }
        
        ResultType res = f_01+f_11;
        if ( res != 0) {
            res = (ResultType) (f_01) / res;
        } else {
            res = 1;
        }
        
        return res;
     };
};

#endif
