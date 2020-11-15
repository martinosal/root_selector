#include <math.h>


int overlap(std::vector<int> x1, std::vector<int> x2){
  int n=0;
  for(unsigned i=0;i<x1.size();i++){
    int a=x1.at(i);
    for(unsigned j=0;j<x2.size();j++){
      if(a==x2.at(j)){
        n++;
      }
    }
  }
  return n;
}


double shrinking_cone_DR(double pT, double m_p1, double m_p2, double m_p3){//pT in GeV !

  double y = m_p1+exp(m_p2 + m_p3 * pT);
  return y;
}
