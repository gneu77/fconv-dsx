#ifndef __STAT_GN
#define __STAT_GN

#include<stdlib.h>
#include<iostream>
#include<vector>
#include<map>
#include<set>
#include<cmath>


using namespace std;


namespace stat_GN {
class VAL_ID {
    public:
        float value;
        int index;
        float rank;
        VAL_ID(float const& v,int const& i):value(v),index(i),rank(999999999.) {}
        ~VAL_ID() {}
        bool equal_val(VAL_ID const& rechts,float const& tolerance) {
            if (fabs(value - rechts.value) <= tolerance) return true;
            else return false;
    }
};


struct VAL_ID_COMP {
  bool operator()(VAL_ID* links,VAL_ID* rechts) const {
      return (links->value < rechts->value);
  }
};


float get_pearson(vector<float>& x_vals,vector<float>& y_vals) {
    if (x_vals.size() != y_vals.size()) {
        cerr << "error:  call to stat_GN::get_pearson with different number of variables (" << x_vals.size() << " x and "
             << y_vals.size() << "y)" << endl;
        exit(1);
    }
    if (x_vals.size() < 2) {
        cerr << "error:  call to stat_GN::get_pearson with less than 2 variables" << endl;
        exit(1);
    }

    float x_mean = 0.;
    float y_mean = 0.;
    for (unsigned int i=0; i<x_vals.size(); ++i) {
        x_mean += x_vals[i];
        y_mean += y_vals[i];
    }
    x_mean /= x_vals.size();
    y_mean /= y_vals.size();
    float sum1 = 0.;
    float sum2 = 0.;
    float sum3 = 0.;
    for (unsigned int i=0; i<x_vals.size(); ++i) {
        sum1 += (x_vals[i]-x_mean) * (y_vals[i]-y_mean);
        sum2 += (x_vals[i]-x_mean) * (x_vals[i]-x_mean);
        sum3 += (y_vals[i]-y_mean) * (y_vals[i]-y_mean);
    }
    sum2 *= sum3;
    return (sum1 / sqrt(sum2));
}


float get_spearman(vector<float>& x_vals,vector<float>& y_vals,float x_tol = 0.000001,float y_tol = 0.000001) {
    //! x_tol * (x_max - x_min) = Toleranzwert = Maximal so gross darf die Differenz zwischen 2 Werten sein,
    //!                                          damit diese als gleich betrachtet werden (=> gleicher Rang)
    if (x_vals.size() != y_vals.size()) {
        cerr << "error:  call to stat_GN::get_spearman with different number of variables (" << x_vals.size() << " x and "
             << y_vals.size() << "y)" << endl;
        exit(1);
    }
    if (x_vals.size() < 2) {
        cerr << "error:  call to stat_GN::get_spearman with less than 2 variables" << endl;
        exit(1);
    }

    multiset<VAL_ID*,VAL_ID_COMP> x_sort;
    multiset<VAL_ID*,VAL_ID_COMP> y_sort;
    float x_min = x_vals[0];
    float x_max = x_vals[0];
    float y_min = y_vals[0];
    float y_max = y_vals[0];
    for (unsigned int i=0; i<x_vals.size(); ++i) {
        x_sort.insert(new VAL_ID(x_vals[i],i));
        if (x_vals[i] < x_min) x_min = x_vals[i];
        if (x_vals[i] > x_max) x_max = x_vals[i];
        y_sort.insert(new VAL_ID(y_vals[i],i));
        if (y_vals[i] < y_min) y_min = y_vals[i];
        if (y_vals[i] > y_max) y_max = y_vals[i];
    }
    float x_range = x_max - x_min;
    float y_range = y_max - y_min;
    x_tol *= x_range;
    y_tol *= y_range;

    int xr = 1;
    vector<VAL_ID*> xv;
    for (multiset<VAL_ID*>::iterator it=x_sort.begin(); it!=x_sort.end(); ++it) {
        (*it)->rank = xr;
        xv.push_back(*it);
        ++xr;
    }
    int yr = 1;
    vector<VAL_ID*> yv;
    for (multiset<VAL_ID*>::iterator it=y_sort.begin(); it!=y_sort.end(); ++it) {
        (*it)->rank = yr;
        yv.push_back(*it);
        ++yr;
    }

    unsigned int last_sim = 0;
    unsigned int hsim = 0;
    for (unsigned int i=1; i<xv.size(); ++i) {
        if (xv[i]->equal_val(*(xv[last_sim]),x_tol)) {
            hsim = i;
        } else {
            if (hsim > last_sim) {
                float mr = 0.;
                for (unsigned int j=last_sim; j<= hsim; ++j) mr += xv[j]->rank;
                mr /= 1 + hsim - last_sim;
                for (unsigned int j=last_sim; j<= hsim; ++j) xv[j]->rank = mr;
            }
            last_sim = i;
            hsim = i;
        }
    }
    if (hsim > last_sim) {
        float mr = 0.;
        for (unsigned int j=last_sim; j<= hsim; ++j) mr += xv[j]->rank;
        mr /= 1 + hsim - last_sim;
        for (unsigned int j=last_sim; j<= hsim; ++j) xv[j]->rank = mr;
    }

    last_sim = 0;
    hsim = 0;
    for (unsigned int i=1; i<yv.size(); ++i) {
        if (yv[i]->equal_val(*(yv[last_sim]),y_tol)) {
            hsim = i;
        } else {
            if (hsim > last_sim) {
                float mr = 0.;
                for (unsigned int j=last_sim; j<= hsim; ++j) mr += yv[j]->rank;
                mr /= 1 + hsim - last_sim;
                for (unsigned int j=last_sim; j<= hsim; ++j) yv[j]->rank = mr;
            }
            last_sim = i;
            hsim = i;
        }
    }
    if (hsim > last_sim) {
        float mr = 0.;
        for (unsigned int j=last_sim; j<= hsim; ++j) mr += yv[j]->rank;
        mr /= 1 + hsim - last_sim;
        for (unsigned int j=last_sim; j<= hsim; ++j) yv[j]->rank = mr;
    }

    map<int,VAL_ID*> xw;
    map<int,VAL_ID*> yw;
    for (unsigned int i=0; i<xv.size(); ++i) {
        xw[xv[i]->index] = xv[i];
        yw[yv[i]->index] = yv[i];
    }

    vector<float> xp;
    vector<float> yp;
    for (map<int,VAL_ID*>::iterator it=xw.begin(); it!=xw.end(); ++it) {
        xp.push_back(it->second->rank);
    }
    for (map<int,VAL_ID*>::iterator it=yw.begin(); it!=yw.end(); ++it) {
        yp.push_back(it->second->rank);
    }

    for (unsigned int i=0; i<xv.size(); ++i) {
        delete xv[i];
        delete yv[i];
    }

    return (get_pearson(xp,yp));
}
}

#endif

