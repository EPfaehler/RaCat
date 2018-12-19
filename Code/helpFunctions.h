#ifndef HELPFUNCTIONS_H_INCLUDED
#define HELPFUNCTIONS_H_INCLUDED

#include <iostream>
using namespace std;



template <typename T>
class square_accumulate
{
public:
    square_accumulate(void) :
      _sum(0)
      {
      }

      const T& result(void) const
      {
          return _sum;
      }

      void operator()(const T& val)
      {
          _sum += val * val;
      }

private:
    T _sum;
};

template <typename T>
class sum_absol_value
{
public:
    sum_absol_value(void) :
      _sum(0)
      {
      }

      const T& result(void) const
      {
          return _sum;
      }

      void operator()(const T& val)
      {
          if(val>0){
          _sum += val;
          }
          else{
          _sum -= val;
          }
      }

private:
    T _sum;
};

template <typename T>
class sum_robust
{
public:
    sum_robust(void) :
      _sum(0)
      {
      }

      const T& result(void) const
      {
          return _sum;
      }

      void operator()(const T& val)
      {
          if(val<0){
          _sum += val;
          }
          else{
          _sum -= val;
          }
      }

private:
    T _sum;
};

#endif //HELPFUNCTIONS_H_INCLUDED




