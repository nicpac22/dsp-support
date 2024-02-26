/* 
 * Copyright (c) 2013 Nick Xenias
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU Lesser General Public License as   
 * published by the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * Lesser General Lesser Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
 
 // Unit tests for constel.h
 // compile with:
 //   g++ -o constel_test constel_test.cc -I../include -O3
 //   icpc -o constel_test constel_test.cc -I../include -O3

#include <vector>
#include "constel.h"

using std::cerr;

// Given 2 vectors, compares their values for equality and returns the index
// of the first element that does not match or -1 if all elements match
template <typename T>
void compare(vector<T> &in1, vector<T> &in2)
{
  int retval = -1;
  if(in1.size() != in2.size()) return min(in1.size(), in2.size());
  for(int i=0; i<in1.size(); ++i)
  {
    if(in1[i] != in2[i])
    {
      retval = i;
      break;
    }
  }
  return retval;
}

void main()
{
  int retval = 0;
  int cval = 0;
  // BPSK tests
  vector<char> ibits(2) = {0,1};
  vector<int> oidx(2);
  vector<int> cidx(2) = {0,1};
  retval

  // QPSK tests

  // 8-PSK tests

  // 16-QAM tests

  // 64-QAM tests
}