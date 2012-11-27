#ifndef CANDIDATE_H
#define CANDIDATE_H

#include "TObject.h"
#include <vector>

using namespace std;

class Candidate : public TObject {
public:
 float chi2, mt2w, mt2bl, mt2b;
 int j1, j2, bi, oi;
 float k1, k2;
 bool match;

 ClassDef(Candidate, 2)
};

typedef vector<Candidate> CANDIDATES;

#endif

