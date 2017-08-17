#include <vector>

#include "basicFunctions.h"

double BF_order1(unsigned int dim, unsigned int num, double u, double v)
{
  double value = 0.;
  switch (dim) {
    case 0:
      /*
       *  phi_0 = 1;
       */
      switch (num) {
        case 0:
          value = 1.;
          break;
        default:
          break;
      }
      break;
    case 1:
      /*
       *  phi_0 = (1-u)/2;
       *  phi_1 = (1+u)/2;
       *
       *  (0)                 (1)
       *  -|---------|---------|-->u
       *   -1        0         1
       */
      switch (num) {
        case 0:
          value = (1-u)/2;
          break;
        case 1:
          value = (1+u)/2;
          break;
        default:
          break;
      }
      break;
    case 2:
      /*
       *  phi_0 = 1-u-v;
       *  phi_1 = u;
       *  phi_2 = v;
       *
       * v ^
       * 1 |- (2)
       *   |  \__
       *   |     \__
       *   |        \__
       *   |           \__
       *   |              \__
       *   |(0)              \ (1)
       *  -|-----------------|-->u
       *   0                 1
       */
      switch (num) {
        case 0:
          value = 1-u-v;
          break;
        case 1:
          value = u;
          break;
        case 2:
          value = v;
          break;
        default:
          break;
      }
      break;
    default:
      break;
  }
  
  return value;
}

std::vector<double> BFgrad_order1(unsigned int dim, unsigned int num, double u, double v)
{
  std::vector<double> value(2, 0.);
  switch (dim) {
    case 0:
      /*
       *  phi_0 = 0;
       */
      switch (num) {
        case 0:
          value[0] = 0.;
          break;
        default:
          break;
      }
      break;
    case 1:
      /*
       *  phi_0 = (1-u)/2;
       *  phi_1 = (1+u)/2;
       *
       *  -|---------|---------|-->u
       *   -1        0         1
       */
      switch (num) {
        case 0:
          value[0] = -1./2;
          break;
        case 1:
          value[0] = 1./2;
          break;
        default:
          break;
      }
      break;
    case 2:
      /*
       *  phi_0 = 1-u-v;
       *  phi_1 = u;
       *  phi_2 = v;
       *
       * v ^
       * 1 |- (2)
       *   |  \__
       *   |     \__
       *   |        \__
       *   |           \__
       *   |              \__
       *   |(0)              \ (1)
       *  -|-----------------|-->u
       *   0                 1
       */
      switch (num) {
        case 0:
          value[0] = -1.;
          value[1] = -1.;
          break;
        case 1:
          value[0] = 1.;
          value[1] = 0.;
          break;
        case 2:
          value[0] = 0.;
          value[1] = 1.;
          break;
        default:
          break;
      }
      break;
    default:
      break;
  }
  
  return value;
}
