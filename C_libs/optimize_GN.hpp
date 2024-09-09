
//============================================================================
// optimize_GN.hpp -*- C++ -*-; optimization and curve fitting
//
// Copyright (C) 2010, 2011 Gerd Neudert
//
// Generic implementation for particle swarm optimization (PSO) and the
// downhill-simplex algorithm by Nelder and Mead.
// Generic curve fitting is realized using PSO. Variance estimates are
// obtained by simple bootstrapping.
// Please read the examples below and have a look at the declaration
// of the constructors and the comments there.
//
// All classes and functions are templates and can be used with
// types float or double (or other floating point types).
//
//      class PS_OPT<T,U>      for particle swarm optimization
//      class NM_OPT<T,U>      for Nelder Mead optimization
//      class CURVE_FIT<T,U>   for curve fitting (using PS_OPT)
//      class PARTICLE<T>      represents a point in search space
//      namespace opt_gn_ns    some little helpers
//
//-------------------------
// CURVE_FIT example usage:
//
//  -> Let us say We have 3 datapoints as dependant variables of an
//     independant variable x: f(x=0)=2, f(x=2)=3 and f(x=5)=4.5
//     Now We like to fit a line through these points, hence We like
//     to find parameters 'a' and 'b' that minimize  f(x) = ax + b
//     The following code will do so and afterwards I will explain it:
//
//  1 int const n_vars = 1; int const n_facs = 2; int const n_data = 3;
//  2 float min[n_facs] = {-5000.,-5000.}; float max[n_facs] = {5000.,5000.};
//  3 float data[n_data] = {2.,3.,4.5};
//  4 float** vars = new float*[n_data];
//  5 for (int i=0; i<n_data; ++i) vars[i] = new float[n_vars];
//  6 vars[0][0] = 0.; vars[1][0] = 2.; vars[2][0] = 5.;
//  7 opt_gn_ns::float_model_func pmf = opt_gn_ns::poly<float,n_vars,1>;
//  8 CURVE_FIT<float,opt_gn_ns::float_model_func> myfit(n_vars,n_facs,n_data,
//                                                       data,vars,pmf,0,
//                                                       min,max);
//  9 myfit.fit(opt_gn_ns::high);
// 10 cout << "f(x) = " << myfit.get_factors()[0] << "x + "
//         << myfit.get_factors()[1] << endl;
//
//  -> In the first line We set some variables. We have only one independant
//     variable (x), We search two factors (a,b) and we have three datapoints.
//  -> In line 2 We limit the search space for the factors.
//  -> In line 3 We store the three known datapoints.
//  -> In line 4 We create a pointer to the variables for each datapoint.
//  -> In line 5 We reserve space for those variables.
//  -> In line 6 We store those variables.
//  -> In line 7 We make a function pointer to our model function. In the
//     namespace opt_gn_ns there are typedefs for common function pointers
//     and there are also generic definitions for some model functions.
//     opt_gn_ns::poly<T,int,int> is a function template that covers many
//     polynomial functions. The first int specifies the number of independant
//     variables and the second the maximum exponent (So for e.g. for the
//     function f(x,y) = ax^3 + bx^2 + cx + dy^3 + ey^2 + fy + g one could
//     use opt_gn_ns::poly<float,2,3>).
//  -> In line 8 We instanciate a CURVE_FIT object. The second template
//     parameter specifies the function pointer type. This is necessary,
//     because it is also possible to use functors instead of function
//     pointers, which makes the usage much more generic.
//  -> In line 9 We apply the fitting with high quality (which is the default,
//     but you may choose between low, medium, high and very_high).
//  -> In line 10 We finally print our result, which hopefully is f(x)=0.5x+2
//
//----------------------
// PS_OPT example usage:
//  -> Though the next example is much shorter, in the case of PSO it is
//     much more important to read the documentation in the constructor
//     declaration very careful.
//  -> The following example minimizes the rosenbrock banana function,
//     which is also predefined in the opt_gn_ns namespace:
//
//  1 float min[2] = {-5000.,-5000.}; float max[2] = {5000.,5000.};
//  2 opt_gn_ns::float_fit_func ff = opt_gn_ns::rosenbrock<float>;
//  3 PS_OPT<float,opt_gn_ns::float_fit_func> myps(2,ff,0,min,max,0,0,500);
//  4 myps.minimize(3);
//  5 cout << "f(" << myps.get_best_position()[0] << ","
//         << myps.get_best_position()[1]
//         << ") = " << myps.get_best_fitness() << endl;
//
//  -> In line 1 We limit the search space for both variables.
//  -> In line 2 We create a function pointer to the rosenbrock function.
//  -> In line 3 We instanciate a PS_OPT object. For the second template
//     parameter, please the the CURVE_FIT example.
//  -> In line 4 We apply the minimization and We specify to make continuous
//     global minimizations until 3 subsequent runs yield no improvement.
//  -> In line 5 We finally print our result, which hopefully is f(1,1) = 0
//============================================================================

#ifndef __OPTIMIZE_GN
#define __OPTIMIZE_GN

#include<stdint.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<cmath>
#include<iostream>
#include<vector>


using namespace std;


namespace opt_gn_ns {
    // Konstanten:
    float const max_fitness = 999999999999999.;
    float const min_fitness = -999999999999999.;
    
    int const low = 0;
    int const medium = 1;
    int const high = 2;
    int const very_high = 3;

    float const sqrt2pi = sqrt(2.*acos(-1.));

    // Variablen:
    unsigned int n_dimensions = 0;

    // Typedefs:
    typedef float (*float_fit_func)(float const *const); // Funktionszeigertyp fuer PS_OPT und NM_OPT
    typedef double (*double_fit_func)(double const *const);
    typedef float (*float_model_func)(float const *const,float const *const); // Funktionszeigertyp fuer CURVE_FIT
    typedef double (*double_model_func)(double const *const,double const *const);
    
    // nuetzliche Funktionen:
    template<class T>
    T square_dist(T const *const p1,T const *const p2) {
        T dst = 0.;
        for (unsigned int i=0; i<opt_gn_ns::n_dimensions; ++i) {
            dst += (p1[i]-p2[i]) * (p1[i]-p2[i]);
        }
        return dst;
    }
    template<class T>
    T euclid_dist(T const *const p1,T const *const p2) {
        // ungewichtete euklidische Distanz als einfaches Abstandsmass zwischen Partikeln
        return sqrt(opt_gn_ns::square_dist<T>(p1,p2));
    }

    // Einfache Fitness-Funktionen zum Testen der Optimierungsverfahren:
    template<class T>
    T rosenbrock(T const *const vals) {
        // Zweidimensionale Variante von Rosenbrocks Bananenfunktion als Test
        // fuer die Minimierung:
        // Minimum bei f(1,1) = 0   -> Von diesem Punkt ausgehend ein
        // lang gezogenes Tal.
        return (100.*(vals[1]-vals[0]*vals[0])*(vals[1]-vals[0]*vals[0])
                +(1-vals[0])*(1-vals[0]));
    }
    template<class T>
    T himmelblau(T const *const vals) {
        // Zweidimensionale Himmelblau-Funktion als Test fuer Minimierung:
        // 4 globale Minima:
        // f(3.0,2.0) = 0.0
        // f(-2.805118, 3.131312) = 0.0
        // f(-3.779310, -3.283186) = 0.0
        // f(3.584428, -1.848126) = 0.0
        // ein lokales Maximum:
        // f(-0.270844,-0.923038) = 181.616 \quad
        T t1 = vals[0]*vals[0]+vals[1]-11;
        t1 *= t1;
        T t2 = vals[0]+vals[1]*vals[1]-7;
        t2 *= t2;
        return (t1+t2);
    }
    template<class T>
    T helical_valley(T const *const vals) {
        // 3D-Funktion mit Minimum bei:
        // f(1.,0.,0.) = 0.
        T t1 = opt_gn_ns::max_fitness;
        if (vals[0] != 0.) t1 = atan(vals[1]/vals[0]);
        if (vals[0] < 0.) t1 += acos(-1.);
        t1 /= 2.*acos(-1.);
        T t2 = sqrt(vals[0]*vals[0] + vals[1]*vals[1]) - 1.; t2 *= t2;
        T t3 = vals[2] - 10.*t1; t3 *= t3;
        return (100.*(t2+t3) + vals[2]*vals[2]);
    }

    // Einfache Funktionen zum Kurven fitten:
    template<class T,int dim,int max_exp>
    T poly(T const *const x,T const *const f) {
        // Beispiel: dim = 2  und  max_exp = 3
        // => f(x_0,x_1) = f_0 * x_0^3 + f_1 * x_0^2 + f_2 * x_0
        //               + f_3 * x_1^3 + f_4 * x_1^2 + f_5 * x_1
        //               + f_6
        // => Es muss also immer  dim  Variablen x  und  dim*max_exp+1  Faktoren f geben
        // - Kennt man die Funktion genauer, sollte natuerlich lieber eine spezialisierte
        //   Version verwendet werden, um Faktoren zu sparen (die sich zu Null ergeben sollten).
        T res = 0.; int nf = 0;
        for (int i=0; i<dim; ++i) {
            T xx = 1.;
            for (int j=0; j<max_exp; ++j) {
                xx *= x[i];
                res += f[nf]*xx;
                ++nf;
            }
        }
        res += f[nf];
        return res;
    }
    template<class T,int dim>
    T gauss(T const *const x,T const *const f) {
        // f[0] = sigma_1, f[1] = mu_1, ...
        // Beispiel: dim = 2
        // => f(x_0,x_1) = 1/(f_0 * sqrt(2PI)) * exp(-0.5*((x_0-f_1)/f_0)^2)
        //               + 1/(f_2 * sqrt(2PI)) * exp(-0.5*((x_1-f_3)/f_2)^2)
        //               + f_4
        T res = 0.;
        for (int i=0; i<dim; ++i) {
            T tq = (x[i]-f[i+1]) / f[i];
            tq *= tq; tq *= -0.5;
            tq = exp(tq);
            tq /= f[i] * sqrt2pi;
            res += tq;
        }
        res += f[2*dim];
        return res;
    }
    template<class T>
    T logistic(T const *const x,T const *const f) {
        // Logistische Funktion: f(x) = [f_0 / (1 + exp(f_1*(f_2-x)))] + f_3
        // => obere Grenze der Sigmoide     = f_0 + f_3
        //    Steigung im Wendepunkt        = f_0 * f_1 / 4.
        //    Wendepunkt   bei            x = f_2
        T res = f[0];
        res /= 1 + exp(f[1]*(f[2]-x[0]));
        res += f[3];
        return res;
    }
    template<class T>
    T dose_response(T const *const x,T const *const f) {
        // Logistische Funktion: f(x) = [f_0 / (1 + exp(f_1*log(x/f_2)))] + f_3
        // => Range      = f_0
        //    Background = f_3
        //    IC50       = f_2
        //    Slope      = f_1
        if (fabs(x[0]) < 0.000001) return (f[0]+f[3]);
        else if (fabs(f[2]) < 0.000001) return f[3];        
        T res = f[0];
        T help = x[0];
        help /= f[2];
        res /= 1 + exp(f[1]*log(help));
        res += f[3];
        return res;
    }
    template<class T>
    T michaelis_menten(T const *const x,T const *const f) {
        // Michaelis Menten Funktion: f(x) = f_0 * x / (f_1 + x)
        // => v_max = f_0
        //    Km    = f_1
        T res = f[0];
        res *= x[0];
        res /= f[1] + x[0];
        return res;
    }
}


template<class T> class PARTICLE;
template<class T> ostream& operator<<(ostream &os,PARTICLE<T> const& ref);


template<class T>
class PARTICLE {
    // Klasse zur Repraesentation eines Individuums im n-dimensionalen Suchraum:
    public:
        int pos_size;
        T fitness;
        T best_fitness;
        T* position;
        T* best_position;
        T* speed;
        int id;
        static uint64_t global_id;
        PARTICLE(int const& p_size):pos_size(p_size),fitness(opt_gn_ns::max_fitness),
                 best_fitness(opt_gn_ns::max_fitness) {
            position = new T[pos_size];
            memset(position,0,pos_size*sizeof(T));
            best_position = new T[pos_size];
            speed = new T[pos_size];
            memset(speed,0,pos_size*sizeof(T));
            id = global_id; ++global_id;
        }
        PARTICLE(int const& p_size,T const *const ref_pos,T const *const ref_speed=0):pos_size(p_size),
                 fitness(opt_gn_ns::max_fitness),best_fitness(opt_gn_ns::max_fitness) {
            position = new T[pos_size];
            memcpy(position,ref_pos,pos_size*sizeof(T));
            best_position = new T[pos_size];
            speed = new T[pos_size];
            if (ref_speed) memcpy(speed,ref_speed,pos_size*sizeof(T));
            else memset(speed,0,pos_size*sizeof(T));
            id = global_id; ++global_id;
        }
        PARTICLE(PARTICLE<T> const& p):pos_size(p.pos_size),fitness(p.fitness),best_fitness(p.best_fitness) {
            position = new T[pos_size];
            memcpy(position,p.position,pos_size*sizeof(T));
            best_position = new T[pos_size];
            memcpy(best_position,p.best_position,pos_size*sizeof(T));
            speed = new T[pos_size];
            memcpy(speed,p.speed,pos_size*sizeof(T));
            id = global_id; ++global_id;
        }
        PARTICLE(PARTICLE<T> const *const p):pos_size(p->pos_size),fitness(p->fitness),best_fitness(p->best_fitness) {
            position = new T[pos_size];
            memcpy(position,p->position,pos_size*sizeof(T));
            best_position = new T[pos_size];
            memcpy(best_position,p->best_position,pos_size*sizeof(T));
            speed = new T[pos_size];
            memcpy(speed,p->speed,pos_size*sizeof(T));
            id = global_id; ++global_id;
        }
        PARTICLE<T>& operator=(PARTICLE<T> const& p) {
            pos_size = p.pos_size;
            fitness = p.fitness;
            best_fitness = p.best_fitness;
            position = new T[pos_size];
            memcpy(position,p.position,pos_size*sizeof(T));
            best_position = new T[pos_size];
            memcpy(best_position,p.best_position,pos_size*sizeof(T));
            speed = new T[pos_size];
            memcpy(speed,p.speed,pos_size*sizeof(T));
            id = global_id; ++global_id;
            return *this;
        }
        ~PARTICLE() {delete[] position; delete[] best_position; delete[] speed;}
        friend ostream &operator<< <>(ostream& os,PARTICLE<T> const& ref);
};
template<class T> uint64_t PARTICLE<T>::global_id = 0;


template<class T>
ostream& operator<<(ostream &os,PARTICLE<T> const& ref) {
    os << "(" << ref.position[0];
    for (int i=1; i<ref.pos_size; ++i) os << "," << ref.position[i];
    os << ") => " << ref.fitness;
    return os;
}


template<class T,class U>
class PS_OPT {
    // Klasse fuer Partikelschwarmoptimierung:
    // - Die meisten Parameter koennen im Konstruktor gesetzt werden,
    //   aber einige sind nur ueber set-Methoden zugaenglich.
    // - T sollte ein Typ fuer Fliesskommazahlen sein.
    // - U kann der Typ fuer einen Funktionspointer auf die Fitnessfunktion
    //   oder aber ein entsprechender Funktor sein.
    private:
        unsigned int n_dimensions;
        vector<PARTICLE<T>* > particles;
        vector<PARTICLE<T>* > trajectory;
        PARTICLE<T>* best_particle;
        PARTICLE<T>* total_best_particle;

        U& fitness_function; // Funktionszeiger oder Funktor
        T* min_val;
        T* max_val;
        int* is_cyclic;
        T* max_speed;
        unsigned int n_particles;
        unsigned int max_iterations;
        unsigned int max_no_improve;
        T (*distance_function)(T const *const,T const *const);
        T max_neighbor_dist;

        //----------------------------------------------------------------------------
        // Parameter die nicht ueber Konstruktor, sondern set-Methoden gesetzt werden:
        T min_percent_change; // Aenderungen < (min_percent_change * fitness) werden nicht als Verbesserung gewertet
        T initial_weight;     // Mit diesem Faktor geht die momentane Geschwindigkeit eines Partikels in die neue
                              //  Geschwindigkeit mit ein. (also eine Art Traegheitsfaktor)
        T weight_decrease;    // Mit diesem Faktor aendert sich das initial_weight pro Iteration.
        T c1_par;             // Mit diesem Faktor geht der Gradient zum lokal besten Punkt in die neue Geschw. ein.
        T c2_par;             // Mit diesem Faktor geht der Gradient zum global besten Punkt in die neue Geschw. ein.
        T c3_par;             // Mit diesem Faktor geht der Gradient zum besten Punkt der Nachbarschaft in die neue Geschw. ein.
        bool track_best;      // Kann auf true gesetzt werden, um eine best_particle-Trajektorie zu speichern
        //----------------------------------------------------------------------------

        int* neighbor_matrix;
        bool search_min;
        int n_start_pos;

        T optimize(int const& max_no_glob_improve=0);
        void generate_particles();
        void initialize_speed();
        void update_neighbor_matrix();
    public:
        PS_OPT(unsigned int const& n_dim,         // Dimensionen

               U& ff,                             // Zeiger auf Fitnessfunktion oder Funktor. (Siehe CURVE_FIT zur Implementierung
                                                  // eines Funktors, der mit nicht-statischen Membern arbeitet.)

               T const *const start_pos=0,        // Eventuelle Startposition (Es koennen weitere Startpositionen mit add_start_position()
                                                  // zugefuegt werden. Aus den vorgegebenen Startpositionen werden initiale Partikel
                                                  // erzeugt. Vor dem Optimieren werden weitere Partikel per Zufall erzeugt, bis n_part
                                                  // Partikel vorhanden sind. Die Startgeschwindigkeiten werden fuer alle Partikel
                                                  // randomisiert, also auch die vorgegebenen.)

               T const *const min=0,              // Minimaler erlaubter Wert fuer jede Dimension.
               T const *const max=0,              // Maximaler erlaubter Wert (min und max begrenzen also den Suchraum. Wird hier nichts
                                                  // vorgegeben ist der Default -1E9 bis 1E9 fuer jede Dimension.)

               int const *const is_cyc=0,         // Cyclische Werte? (z.B. 0 bis 2PI mit 2PI=0  => Wird hier fuer eine Dimension
                                                  // Null angegeben ist sie nicht cyclisch. Bei einem Wert ungleich Null ist sie
                                                  // cyclisch, d.h. min und max entsprechen der gleichen Position innerhalb dieser
                                                  // Dimension. Mit min=0 und max=7 wäre z.B. die kuerzeste Distanz zwischen den
                                                  // Punkten 1 und 5 nicht 4, sondern -3)

               T const *const max_sp=0,           // Maximal erlaubte Geschwindigkeit (Spezifiziert wie weit sich ein Partikel maximal
                                                  // pro Iteration innerhalb der jeweiligen Dimension bewegen darf. Wenn hier kein Limit
                                                  // vorgegeben wird, dann wird per Default der mittlere interpartikulaere Abstand
                                                  // innerhalb der jeweiligen Dimension genommen.)

               unsigned int const& n_part=0,      // Wieviele Partikel sollen generiert werden? Bei einem Wert von Null werden
                                                  // per Default  n_dim * 30  Partikel genommen.

               unsigned int const& max_iter=500,  // Wieviele Iterationen maximal (Nach dieser Anzahl Iterationen wird abgebrochen,
                                                  // wenn vorher kein anderes Abbruchkriterium erfuellt wird.)
               unsigned int const& max_no_imp=25, // Vorzeitiger Abbruch, wenn entsprechend viele Iterationen keine Verbesserung
                                                  // bringen (Verbesserungen < 0.1 Prozent werden nicht als Verbesserungen gezaehlt.
                                                  // Dieser Prozentwert kann ueber set_min_percent_change() geaendert werden.)
                
               T (*df)(T const *const,            // Zeiger auf Distanzfunktion (Wenn neben global und lokal bester Position auch die
                       T const *const)=           // beste Position benachbarter Partikel beruecksichtigt werden soll (max_nd>0), dann
                       opt_gn_ns::euclid_dist<T>, // wird eine Distanzfunktion benoetigt um den Abstand zwischen 2 Partikeln zu bestimmen.
                                                  // Default ist die einfache ungewichtete euklidische Distanz.)
               T const& max_nd=-1.);              // Maximaler Abstand fuer Nachbarschaft (Ein Wert <= 0. bedeutet, dass keine Nachbarschaft
                                                  // beruecksichtigt wird. Fuer Werte > 0. werden alle Partikel mit einem Abstand <= max_nd
                                                  // als Nachbarn betrachtet. ACHTUNG: Dies vermeidet schnelle Konvergenz gegen lokale Minima
                                                  // und verbessert somit die Ergebnisse fuer Energiehyperflaechen mit sehr vielen lokalen
                                                  // Minima. Der Rechenaufwand fuer eine Iteration steigt jedoch erheblich, da die
                                                  // interpartikulaeren Abstaende in jeder Iteration neu bestimmt werden muessen!)

        ~PS_OPT();

        void add_start_position(T const *const start_pos); // zusaetzliche Startpositionen vorgeben
        void set_min_percent_change(T const& val);
        void set_initial_weight(T const& val);
        void set_weight_decrease(T const& val);
        void set_c1_par(T const& val);
        void set_c2_par(T const& val);
        void set_c3_par(T const& val);
        void set_track_best(bool const& val);

        T minimize(int const& max_no_glob_improve=0); // Wird hier ein Wert>0 mitgegeben, so wird die globale Suche mehrfach gestartet.
        T maximize(int const& max_no_glob_improve=0); // Ein Wert von 3 bedeutet z.B., dass solange die Minimierung neu gestartet wird,
                                                      // bis 3 aufeinander folgende Minimierungen keine Verbesserung bringen.

        T const *const get_best_position();
        T const get_best_fitness();
        vector<PARTICLE<T>* >& get_trajectory();
};


template<class T,class U>
class NM_OPT {
    // Klasse fuer Downhill Simplex Optimierung nach Nelder-Mead:
    // - Die meisten Parameter koennen im Konstruktor gesetzt werden,
    //   aber einige sind nur ueber set-Methoden zugaenglich.
    // - T sollte ein Typ fuer Fliesskommazahlen sein.
    // - U kann der Typ fuer einen Funktionspointer auf die Fitnessfunktion
    //   oder aber ein entsprechender Funktor sein.
    private:
        unsigned int n_dimensions;
        unsigned int n_points;
        vector<PARTICLE<T>* > particles;
        vector<PARTICLE<T>* > trajectory;
        PARTICLE<T>* total_best_particle;

        U& fitness_function; // Funktionszeiger oder Funktor
        T* min_val;
        T* max_val;
        int* is_cyclic;
        unsigned int max_iterations;
        T min_square_dist;
        bool use_start_pos;

        //----------------------------------------------------------------------------
        // Parameter die nicht ueber Konstruktor, sondern set-Methoden gesetzt werden:
        T alpha;              // Faktor fuer Reflektion
        T beta;               // Faktor fuer Kontraktion
        T gamma;              // Faktor fuer Expansion
        T delta;              // Faktor fuer Reduktion
        bool track_best;      // Kann auf true gesetzt werden, um eine best_particle-Trajektorie zu speichern
        //----------------------------------------------------------------------------

        bool search_min;

        T optimize();
        void generate_start_positions();
    public:
        NM_OPT(unsigned int const& n_dim,         // Dimensionen
               U& ff,                             // Zeiger auf Fitnessfunktion oder Funktor. (Siehe CURVE_FIT zur Implementierung
                                                  // eines Funktors, der mit nicht-statischen Membern arbeitet.)

               T const *const start_pos=0,        // Startposition. Wird kein Start vorgegeben wird (max-min)/2 verwendet
               bool const& use_start=true,        // true  => start_pos wird Teil des initialen Simplex
                                                  // false => start_pos wird Mittelpunkt des initialen Simplex
               
               T const *const min=0,              // Minimaler erlaubter Wert fuer jede Dimension.
               T const *const max=0,              // Maximaler erlaubter Wert (min und max begrenzen also den Suchraum. Wird hier nichts
                                                  // vorgegeben ist der Default -1E9 bis 1E9 fuer jede Dimension.)

               int const *const is_cyc=0,         // Cyclische Werte? (z.B. 0 bis 2PI mit 2PI=0  => Wird hier fuer eine Dimension
                                                  // Null angegeben ist sie nicht cyclisch. Bei einem Wert ungleich Null ist sie
                                                  // cyclisch, d.h. min und max entsprechen der gleichen Position innerhalb dieser
                                                  // Dimension. Mit min=0 und max=7 wäre z.B. die kuerzeste Distanz zwischen den
                                                  // Punkten 1 und 5 nicht 4, sondern -3)

               unsigned int const& max_iter=5000, // Wieviele Iterationen maximal (Nach dieser Anzahl Iterationen wird abgebrochen,
                                                  // wenn vorher kein anderes Abbruchkriterium erfuellt wird.)
               T const& min_sdst=-1.);            // Vorzeitiger Abbruch, wenn das Quadrat der Distanz zwischen bestem und schlechtestem
                                                  // Punkt des Simplex unter diesen Wert sinkt. Beim Default -1. wird automatisch
                                                  // das Quadrat von einem tausendstel der initialen Distanz zwischen bestem und
                                                  // schlechtesten Punkt genommen.
        ~NM_OPT();

        void set_alpha(T const& val);
        void set_beta(T const& val);
        void set_gamma(T const& val);
        void set_delta(T const& val);
        void set_track_best(bool const& val);

        T minimize();
        T maximize();

        T const *const get_best_position();
        T const get_best_fitness();
        vector<PARTICLE<T>* >& get_trajectory();
};


template<class T,class U>
class CURVE_FIT {
    private:
        unsigned int n_dimensions;
        unsigned int n_factors;
        U& model_function;
        unsigned int n_data_points;
        int used_quality;
        T* data_points;
        T** data_vars;
        T* min_val;
        T* max_val;
        T* factors;
        T* start_vals;
        T* bs_means;        
        T* bs_variances;
        bool bootstrap_calculated;
        class SQUARE_DIFF { // Funktor zur Berechnung der quadratischen Abweichung
            private:        // zwischen Modell und Datenpunkten.
                CURVE_FIT* parent;
            public:
                T operator()(T const *const vals);
                SQUARE_DIFF(CURVE_FIT* cf):parent(cf) {} // Zeiger auf CURVE_FIT Instanz, damit auf non-static
                ~SQUARE_DIFF() {}                        // Member zugegriffen werden kann.
        };
        void bootstrap(unsigned int const& n_simulations = 100);
    public:
        CURVE_FIT(unsigned int const& n_dim,      // Dimensionen (z.B. hat f(x)=ax+b die Dimension 1, wobei aber
                                                  // zwei Faktoren als Ergebnis ermittelt werden)

                  unsigned int const& n_facs,     // Anzahl der zu ermittelnden Faktoren (bei f(x)=ax+b ist n_facs=n_dim+1
                                                  // Bei anderen Funktionen kann das natuerlich ganz anders sein

                  unsigned int const& n_data,     // Wieviele Datenpunkte werden uebergeben?
                  T const *const data,            // Die Datenpunkte (z.B.: f(x,y)=10 und f(x,y)=-13
                  T const *const *const vars,     // Die Variablen zu den Datenpunkten (z.B.: x=3,y=5 und x=1,y=-2)
                                                  // vars[i][j] ist Variable j von Funktionswert i

                  U& mf,                          // Zeiger auf Modellfunktion oder Funktor. Das erste Argument zeigt
                                                  // auf die Variablenwerte und das zweite auf die Faktoren.

                  T const *const start_param=0,   // Eventuelle Startparameter (Faktoren)

                  T const *const min=0,           // Minimaler erlaubter Wert fuer jeden Faktor.
                  T const *const max=0);          // Maximaler erlaubter Wert (min und max begrenzen also den Suchraum. Wird hier
                                                  // nichts vorgegeben ist der Default -1E9 bis 1E9 fuer jede Dimension.)
        ~CURVE_FIT();

        T fit(int const& quality=opt_gn_ns::high);   // Es stehen low, medium, high und very_high zur Auswahl (Hoehere Qualitaet
                                                     // bedingt natuerlich hoehere Laufzeit.)
        T const *const get_factors();
        T const get_error();
        T const *const get_bootstrap_variances(unsigned int const& n_simulations = 100);
        T const *const get_bootstrap_means(unsigned int const& n_simulations = 100);
};


//! #################################################################################################################
//! #################################################################################################################


template<class T,class U>
PS_OPT<T,U>::PS_OPT(unsigned int const& n_dim,U& ff,T const *const start_pos,
                    T const *const min,T const *const max,int const *const is_cyc,
                    T const *const max_sp,unsigned int const& n_part,unsigned int const& max_iter,
                    unsigned int const& max_no_imp,T (*df)(T const *const,T const *const),
                    T const& max_nd):n_dimensions(n_dim),fitness_function(ff) {
    opt_gn_ns::n_dimensions = n_dimensions; // wird nur bei max_nd > 0. gebraucht
    if (n_part > 0) n_particles = n_part;
    else n_particles = n_dimensions * 30;
    max_iterations = max_iter;
    max_no_improve = max_no_imp;
    max_neighbor_dist = max_nd;
    distance_function = df;
    n_start_pos = 0;
    if (max_neighbor_dist > 0.) {
        neighbor_matrix = new int[n_particles*n_particles];
    } else neighbor_matrix = 0;

    if (start_pos) {
        particles.push_back(new PARTICLE<T>(n_dimensions,start_pos));
        ++n_start_pos;
    }

    min_val = new T[n_dimensions];
    if (min) memcpy(min_val,min,n_dimensions*sizeof(T));
    else for (unsigned int i=0; i<n_dimensions; ++i) min_val[i] = -999999999.;

    max_val = new T[n_dimensions];
    if (max) memcpy(max_val,max,n_dimensions*sizeof(T));
    else for (unsigned int i=0; i<n_dimensions; ++i) max_val[i] = 999999999.;

    is_cyclic = new int[n_dimensions];
    if (is_cyc) memcpy(is_cyclic,is_cyc,n_dimensions*sizeof(int));
    else for (unsigned int i=0; i<n_dimensions; ++i) is_cyclic[i] = 0;

    max_speed = new T[n_dimensions];
    if (max_sp) memcpy(max_speed,max_sp,n_dimensions*sizeof(T));
    else for (unsigned int i=0; i<n_dimensions; ++i) {
        T inter_dist = (max_val[i] - min_val[i]) / n_particles;
        max_speed[i] = inter_dist;
    }

    //-----------------------------------------------------------------------------------
    // Defaults fuer die Parameter, die nur ueber set-Methoden geaendert werden koennen:
    min_percent_change = 0.001;
    initial_weight = 0.9;
    weight_decrease = 0.99;
    c1_par = 2.0;
    c2_par = 2.0;
    c3_par = 1.0;
    track_best = false;
    //-----------------------------------------------------------------------------------

    search_min = true;
    best_particle = new PARTICLE<T>(n_dimensions);
    best_particle->id = -1;
    total_best_particle = new PARTICLE<T>(n_dimensions);
    total_best_particle->id = -1;

    srand(time(0));
}


template<class T,class U>
PS_OPT<T,U>::~PS_OPT() {
    for (typename vector<PARTICLE<T>* >::iterator it=particles.begin(); it!=particles.end(); ++it) {
        delete *it;
    }
    for (typename vector<PARTICLE<T>* >::iterator it=trajectory.begin(); it!=trajectory.end(); ++it) {
        delete *it;
    }
    delete best_particle;
    delete total_best_particle;
    delete[] min_val;
    delete[] max_val;
    delete[] is_cyclic;
    delete[] max_speed;
    if (neighbor_matrix) delete[] neighbor_matrix;
}


template<class T,class U>
void PS_OPT<T,U>::add_start_position(T const *const start_pos) {
    particles.push_back(new PARTICLE<T>(n_dimensions,start_pos));
    ++n_start_pos;
}


template<class T,class U> void PS_OPT<T,U>::set_min_percent_change(T const& val) {min_percent_change = val / 100.;}
template<class T,class U> void PS_OPT<T,U>::set_initial_weight(T const& val) {initial_weight = val;}
template<class T,class U> void PS_OPT<T,U>::set_weight_decrease(T const& val) {weight_decrease = val;}
template<class T,class U> void PS_OPT<T,U>::set_c1_par(T const& val) {c1_par = val;}
template<class T,class U> void PS_OPT<T,U>::set_c2_par(T const& val) {c2_par = val;}
template<class T,class U> void PS_OPT<T,U>::set_c3_par(T const& val) {c3_par = val;}
template<class T,class U> void PS_OPT<T,U>::set_track_best(bool const& val) {track_best = val;}


template<class T,class U> T const *const PS_OPT<T,U>::get_best_position() {return total_best_particle->position;}
template<class T,class U> T const PS_OPT<T,U>::get_best_fitness() {return total_best_particle->fitness;}
template<class T,class U> vector<PARTICLE<T>* >& PS_OPT<T,U>::get_trajectory() {return trajectory;}


template<class T,class U>
T PS_OPT<T,U>::minimize(int const& max_no_glob_improve) {
    search_min = true;
    return optimize(max_no_glob_improve);
}


template<class T,class U>
T PS_OPT<T,U>::maximize(int const& max_no_glob_improve) {
    search_min = false;
    return optimize(max_no_glob_improve);
}


template<class T,class U>
T PS_OPT<T,U>::optimize(int const& max_no_glob_improve) {
    if (search_min) total_best_particle->fitness = opt_gn_ns::max_fitness;
    else total_best_particle->fitness = opt_gn_ns::min_fitness;
    int glob_no_change = 0;
    while (true) {
        if (search_min) best_particle->fitness = opt_gn_ns::max_fitness;
        else best_particle->fitness = opt_gn_ns::min_fitness;
        generate_particles();
        initialize_speed();

        unsigned int no_iter_improve = 0;
        T w = initial_weight;
        for (unsigned int i=0; i<max_iterations; ++i) {
            ++no_iter_improve;
            // Fitness bestimmen und local- und global-best updaten:
            for (typename vector<PARTICLE<T>* >::iterator it=particles.begin(); it!=particles.end(); ++it) {
                (*it)->fitness = fitness_function((*it)->position);
                if (search_min) {
                    if ((*it)->fitness < (*it)->best_fitness) {
                        (*it)->best_fitness = (*it)->fitness;
                        memcpy((*it)->best_position,(*it)->position,n_dimensions*sizeof(T));
                    }
                    if ((*it)->fitness < best_particle->fitness) {
                        if ((best_particle->fitness-(*it)->fitness) > fabs(min_percent_change*((*it)->fitness))) no_iter_improve = 0;
                        best_particle->fitness = (*it)->fitness;
                        memcpy(best_particle->position,(*it)->position,n_dimensions*sizeof(T));
                        if (track_best && (best_particle->fitness < total_best_particle->fitness)) {
                            trajectory.push_back(new PARTICLE<T>(best_particle));
                        }
                    }
                } else {
                    if ((*it)->fitness > (*it)->best_fitness) {
                        (*it)->best_fitness = (*it)->fitness;
                        memcpy((*it)->best_position,(*it)->position,n_dimensions*sizeof(T));
                    }
                    if ((*it)->fitness > best_particle->fitness) {
                        if (((*it)->fitness-best_particle->fitness) > fabs(min_percent_change*((*it)->fitness))) no_iter_improve = 0;
                        best_particle->fitness = (*it)->fitness;
                        memcpy(best_particle->position,(*it)->position,n_dimensions*sizeof(T));
                        if (track_best && (best_particle->fitness > total_best_particle->fitness)) {
                            trajectory.push_back(new PARTICLE<T>(best_particle));
                        }
                    }
                }
            }

            if (max_neighbor_dist > 0.) update_neighbor_matrix();
            if (no_iter_improve > max_no_improve) break;

            for (typename vector<PARTICLE<T>* >::iterator it=particles.begin(); it!=particles.end(); ++it) {
                // Wenn die Nachbarschaft der Partikel beruecksichtigt werden soll, dann
                // muss der Nachbar mit der besten Fitness bestimmt werden:
                int best_neighbor = (*it)->id;
                if (max_neighbor_dist > 0.) {
                    for (unsigned int k=0; k<particles.size(); ++k) {
                        if (int(k) == (*it)->id) continue;
                        if (neighbor_matrix[(*it)->id*particles.size()+k]) {
                            if (search_min) {
                                if (particles[k]->fitness < particles[best_neighbor]->fitness) {
                                    best_neighbor = k;
                                }
                            } else {
                                if (particles[k]->fitness > particles[best_neighbor]->fitness) {
                                    best_neighbor = k;
                                }
                            }
                        }
                    }
                }

                // Jetzt Geschwindigkeit und Positionen der Partikel updaten:
                for (unsigned int j=0; j<n_dimensions; ++j) {
                    T r1 = rand() / T(RAND_MAX);
                    T r2 = rand() / T(RAND_MAX);

                    T local_diff = (*it)->best_position[j] - (*it)->position[j];
                    T global_diff = best_particle->position[j] - (*it)->position[j];
                    if (is_cyclic[j]) {
                        if (local_diff < 0.) {
                            T tester = max_val[j]-(*it)->position[j] + (*it)->best_position[j]-min_val[j];
                            if (fabs(tester) < fabs(local_diff)) local_diff = tester;
                        } else {
                            T tester = min_val[j]-(*it)->position[j] + (*it)->best_position[j]-max_val[j];
                            if (fabs(tester) < fabs(local_diff)) local_diff = tester;
                        }
                        if (global_diff < 0.) {
                            T tester = max_val[j]-(*it)->position[j] + best_particle->position[j]-min_val[j];
                            if (fabs(tester) < fabs(global_diff)) global_diff = tester;
                        } else {
                            T tester = min_val[j]-(*it)->position[j] + best_particle->position[j]-max_val[j];
                            if (fabs(tester) < fabs(global_diff)) global_diff = tester;
                        }
                    }

                    (*it)->speed[j] = w * (*it)->speed[j] + r1 * c1_par * local_diff + r2 * c2_par * global_diff;

                    if (max_neighbor_dist > 0. && best_neighbor != (*it)->id) {
                        T neighbor_diff = particles[best_neighbor]->position[j] - (*it)->position[j];
                        if (is_cyclic[j]) {
                            if (neighbor_diff < 0.) {
                                T tester = max_val[j]-(*it)->position[j] + particles[best_neighbor]->position[j]-min_val[j];
                                if (fabs(tester) < fabs(neighbor_diff)) neighbor_diff = tester;
                            } else {
                                T tester = min_val[j]-(*it)->position[j] + particles[best_neighbor]->position[j]-max_val[j];
                                if (fabs(tester) < fabs(neighbor_diff)) neighbor_diff = tester;
                            }
                        }
                        T r3 = rand() / T(RAND_MAX);
                        (*it)->speed[j] += r3 * c3_par * neighbor_diff;
                    }

                    if ((*it)->speed[j] > max_speed[j]) (*it)->speed[j] = max_speed[j];
                    else if ((*it)->speed[j] < -max_speed[j]) (*it)->speed[j] = -max_speed[j];

                    (*it)->position[j] += (*it)->speed[j];
                    if ((*it)->position[j] < min_val[j]) {
                        if (is_cyclic[j]) {
                            (*it)->position[j] = max_val[j] + ((*it)->position[j] - min_val[j]);
                        } else (*it)->position[j] = min_val[j];
                    } else if ((*it)->position[j] > max_val[j]) {
                        if (is_cyclic[j]) {
                            (*it)->position[j] = min_val[j] + ((*it)->position[j] - max_val[j]);
                        } else (*it)->position[j] = max_val[j];
                    }
                }
                w *= weight_decrease;
            }
        }

        if (search_min) {
            if (best_particle->fitness < total_best_particle->fitness) {
                total_best_particle->fitness = best_particle->fitness;
                memcpy(total_best_particle->position,best_particle->position,n_dimensions*sizeof(T));
                glob_no_change = 1;
            }
        } else {
            if (best_particle->fitness > total_best_particle->fitness) {
                total_best_particle->fitness = best_particle->fitness;
                memcpy(total_best_particle->position,best_particle->position,n_dimensions*sizeof(T));
                glob_no_change = 1;
            }
        }
        ++glob_no_change;
        if (glob_no_change > max_no_glob_improve) break;
    }

    return total_best_particle->fitness;
}


template<class T,class U>
void PS_OPT<T,U>::generate_particles() {
    for (unsigned int i=n_start_pos; i<particles.size(); ++i) delete particles[i];
    particles.resize(n_start_pos);
    while (particles.size() < n_particles) {
        PARTICLE<T>* new_part = new PARTICLE<T>(n_dimensions);
        for (unsigned int j=0; j<n_dimensions; ++j) {
            T r = rand() / T(RAND_MAX);
            new_part->position[j] = min_val[j] + (r * (max_val[j] - min_val[j]));
        }
        particles.push_back(new_part);
    }
    if (max_neighbor_dist > 0.) {
        // Wichtig: Damit die Nachbarmatrix richtig indiziert wird, muessen die
        // IDs der Partikel hier nochmal gesetzt werden!
        for (unsigned int i=0; i<particles.size(); ++i) {
            particles[i]->id = i;
        }
    }
}


template<class T,class U>
void PS_OPT<T,U>::initialize_speed() {
    for (typename vector<PARTICLE<T>* >::iterator it=particles.begin(); it!=particles.end(); ++it) {
        for (unsigned int j=0; j<n_dimensions; ++j) {
            T r = 1. - (2. * rand() / T(RAND_MAX));
            (*it)->speed[j] = r * max_speed[j];
        }
        if (search_min) (*it)->best_fitness = opt_gn_ns::max_fitness;
        else (*it)->best_fitness = opt_gn_ns::min_fitness;
    }
}


template<class T,class U>
void PS_OPT<T,U>::update_neighbor_matrix() {
    for (unsigned int i=0; i<particles.size(); ++i) {
        for (unsigned int j=i+1; j<particles.size(); ++j) {
            int idx1 = j*particles.size() + i;
            int idx2 = i*particles.size() + j;
            if (distance_function(particles[i]->position,particles[j]->position) > max_neighbor_dist) {
                neighbor_matrix[idx1] = 0; neighbor_matrix[idx2] = 0;
            } else {
                neighbor_matrix[idx1] = 1; neighbor_matrix[idx2] = 1;
            }
        }
    }
}


//! ============================================================================


template<class T,class U>
NM_OPT<T,U>::NM_OPT(unsigned int const& n_dim,U& ff,T const *const start_pos,bool const& use_start,
                    T const *const min,T const *const max,int const *const is_cyc,unsigned int const& max_iter,
                    T const& min_sdst):n_dimensions(n_dim),fitness_function(ff) {
    opt_gn_ns::n_dimensions = n_dimensions;
    n_points = n_dimensions + 1;
    max_iterations = max_iter;
    min_square_dist = min_sdst;
    use_start_pos = use_start;

    min_val = new T[n_dimensions];
    if (min) memcpy(min_val,min,n_dimensions*sizeof(T));
    else for (unsigned int i=0; i<n_dimensions; ++i) min_val[i] = -999999999.;

    max_val = new T[n_dimensions];
    if (max) memcpy(max_val,max,n_dimensions*sizeof(T));
    else for (unsigned int i=0; i<n_dimensions; ++i) max_val[i] = 999999999.;
    
    if (start_pos) total_best_particle = new PARTICLE<T>(n_dimensions,start_pos);
    else {
        total_best_particle = new PARTICLE<T>(n_dimensions);
        for (unsigned int i=0; i<n_dimensions; ++i) total_best_particle->position[i] = (max_val[i]-min_val[i]) / 2.;
    }

    is_cyclic = new int[n_dimensions];
    if (is_cyc) memcpy(is_cyclic,is_cyc,n_dimensions*sizeof(int));
    else for (unsigned int i=0; i<n_dimensions; ++i) is_cyclic[i] = 0;

    //-----------------------------------------------------------------------------------
    // Defaults fuer die Parameter, die nur ueber set-Methoden geaendert werden koennen:
    alpha = 1.0;
    beta = 0.5;
    gamma = 2.0; // muss immer > alpha sein!
    delta = 0.5;
    track_best = false;
    //-----------------------------------------------------------------------------------

    search_min = true;

    srand(time(0));
}


template<class T,class U>
NM_OPT<T,U>::~NM_OPT() {
    for (typename vector<PARTICLE<T>* >::iterator it=particles.begin(); it!=particles.end(); ++it) {
        delete *it;
    }
    for (typename vector<PARTICLE<T>* >::iterator it=trajectory.begin(); it!=trajectory.end(); ++it) {
        delete *it;
    }
    delete total_best_particle;
    delete[] min_val;
    delete[] max_val;
    delete[] is_cyclic;
}


template<class T,class U> void NM_OPT<T,U>::set_alpha(T const& val) {alpha = val;}
template<class T,class U> void NM_OPT<T,U>::set_beta(T const& val) {beta = val;}
template<class T,class U> void NM_OPT<T,U>::set_gamma(T const& val) {gamma = val;}
template<class T,class U> void NM_OPT<T,U>::set_delta(T const& val) {delta = val;}
template<class T,class U> void NM_OPT<T,U>::set_track_best(bool const& val) {track_best = val;}


template<class T,class U> T const *const NM_OPT<T,U>::get_best_position() {return total_best_particle->position;}
template<class T,class U> T const NM_OPT<T,U>::get_best_fitness() {return total_best_particle->fitness;}
template<class T,class U> vector<PARTICLE<T>* >& NM_OPT<T,U>::get_trajectory() {return trajectory;}


template<class T,class U>
T NM_OPT<T,U>::minimize() {
    search_min = true;
    return optimize();
}


template<class T,class U>
T NM_OPT<T,U>::maximize() {
    search_min = false;
    return optimize();
}


template<class T,class U>
T NM_OPT<T,U>::optimize() {
    generate_start_positions();

    T* center = new T[n_dimensions];
    T* tp = new T[n_dimensions];
    T* tp2 = new T[n_dimensions];

    unsigned int best_i = 0;
    unsigned int worst_i = 0; unsigned int second_worst_i = 0;
    for (unsigned int i=0; i<n_points; ++i) {
        particles[i]->fitness = fitness_function(particles[i]->position);
        if (search_min) {
            if (particles[i]->fitness < particles[best_i]->fitness) {
                best_i = i;
            } else if (particles[i]->fitness > particles[worst_i]->fitness) {
                second_worst_i = worst_i;
                worst_i = i;
            }
        } else {
            if (particles[i]->fitness > particles[best_i]->fitness) {
                best_i = i;
            } else if (particles[i]->fitness < particles[worst_i]->fitness) {
                second_worst_i = worst_i;
                worst_i = i;
            }
        }
    }
    // Wenn gleich der erste Punkt der schlechteste war, dann stimmt second_worst nicht.
    // Ausserdem soll auch bei gleichen Fitnesswerten worst und second_worst nie gleich sein!:
    if (worst_i == second_worst_i) {
        if (worst_i == 0) second_worst_i = 1;
        else second_worst_i = 0;
        for (unsigned int i=0; i<n_points && i!=worst_i; ++i) {
            if (search_min) {
                if (particles[i]->fitness > particles[second_worst_i]->fitness) second_worst_i = i;
            } else {
                if (particles[i]->fitness < particles[second_worst_i]->fitness) second_worst_i = i;
            }
        }
    }

    if (min_square_dist < 0.) {
        min_square_dist = opt_gn_ns::euclid_dist(particles[best_i]->position,particles[worst_i]->position);
        min_square_dist /= 1000.;
        min_square_dist *= min_square_dist;
    }

    for (unsigned int iter=0; iter<max_iterations; ++iter) {
        // Update von total_best und eventuell Trajektorie erweitern:
        if (search_min) {
            if (particles[best_i]->fitness < total_best_particle->fitness) {
                if (track_best) trajectory.push_back(new PARTICLE<T>(particles[best_i]));
                total_best_particle->fitness = particles[best_i]->fitness;
                memcpy(total_best_particle->position,particles[best_i]->position,n_dimensions*sizeof(T));
            }
        } else {
            if (particles[best_i]->fitness > total_best_particle->fitness) {
                if (track_best) trajectory.push_back(new PARTICLE<T>(particles[best_i]));
                total_best_particle->fitness = particles[best_i]->fitness;
                memcpy(total_best_particle->position,particles[best_i]->position,n_dimensions*sizeof(T));
            }
        }

        // Auf vorzeitigen Abbruch pruefen:
        if (opt_gn_ns::square_dist(particles[best_i]->position,particles[worst_i]->position) < min_square_dist) {
            //cout << "DEBUG: early break after " << iter << " iterations" << endl;
            break;
        }

        // worst_i an center reflektieren:
        for (unsigned int i=0; i<n_dimensions; ++i) {
            // center ohne worst_i berechnen:
            center[i] = 0.;
            for (unsigned int j=0; j<n_points && j!=worst_i; ++j) {
                center[i] += particles[j]->position[i];
            }
            center[i] /= n_dimensions;
        }
        // Reflektion:
        for (unsigned int i=0; i<n_dimensions; ++i) {
            tp[i] = center[i] + alpha*(center[i]-particles[worst_i]->position[i]);
            if (tp[i] > max_val[i]) {
                if (is_cyclic[i]) {
                    tp[i] = min_val[i] + (tp[i] - max_val[i]);
                } else tp[i] = max_val[i];
            } else if (tp[i] < min_val[i]) {
                if (is_cyclic[i]) {
                    tp[i] = max_val[i] + (tp[i] - min_val[i]);
                } else tp[i] = min_val[i];
            }
        }
        T tp_fitness = fitness_function(tp);

        if ((search_min && (tp_fitness < particles[best_i]->fitness)) ||
            (!search_min && (tp_fitness > particles[best_i]->fitness))) {
            // Expansion:
            for (unsigned int i=0; i<n_dimensions; ++i) {
                tp2[i] = center[i] + gamma*(center[i]-particles[worst_i]->position[i]);
                if (tp2[i] > max_val[i]) {
                    if (is_cyclic[i]) {
                        tp2[i] = min_val[i] + (tp2[i] - max_val[i]);
                    } else tp2[i] = max_val[i];
                } else if (tp2[i] < min_val[i]) {
                    if (is_cyclic[i]) {
                        tp2[i] = max_val[i] + (tp2[i] - min_val[i]);
                    } else tp2[i] = min_val[i];
                }
            }
            T tp2_fitness = fitness_function(tp2);
            if ((search_min && (tp2_fitness < tp_fitness)) ||
                (!search_min && (tp2_fitness > tp_fitness))) {
                // worst_i durch tp2 ersetzen:
                memcpy(particles[worst_i]->position,tp2,n_dimensions*sizeof(T));
                particles[worst_i]->fitness = tp2_fitness;
            } else {
                // worst_i durch tp ersetzen:
                memcpy(particles[worst_i]->position,tp,n_dimensions*sizeof(T));
                particles[worst_i]->fitness = tp_fitness;
            }
            // Indizes umsetzen:
            best_i = worst_i;
            worst_i = second_worst_i;
            for (unsigned int i=0; i<n_points && i!=worst_i; ++i) {
                if (search_min) {
                    if (particles[i]->fitness > particles[second_worst_i]->fitness) second_worst_i = i;
                } else {
                    if (particles[i]->fitness < particles[second_worst_i]->fitness) second_worst_i = i;
                }
            }
        } else {
            if ((search_min && (tp_fitness < particles[second_worst_i]->fitness)) ||
                (!search_min && (tp_fitness > particles[second_worst_i]->fitness))) {
                // worst_i durch tp ersetzen:
                memcpy(particles[worst_i]->position,tp,n_dimensions*sizeof(T));
                particles[worst_i]->fitness = tp_fitness;
            } else {
                // Kontraktion:
                if ((search_min && (tp_fitness < particles[worst_i]->fitness)) ||
                    (!search_min && (tp_fitness > particles[worst_i]->fitness))) {
                    // fuer Kontraktion worst_i durch tp ersetzen:
                    for (unsigned int i=0; i<n_dimensions; ++i) {
                        tp2[i] = tp[i] + beta*(center[i]-tp[i]);
                        // Keine min/max Pruefung bei Kontraktion (Grenzen koennen nicht ueberschritten werden)
                    }
                } else {
                    for (unsigned int i=0; i<n_dimensions; ++i) {
                        tp2[i] = particles[worst_i]->position[i] + beta*(center[i]-particles[worst_i]->position[i]);
                    }
                }
                T tp2_fitness = fitness_function(tp2);
                if ((search_min && (tp2_fitness < particles[worst_i]->fitness)) ||
                    (!search_min && (tp2_fitness > particles[worst_i]->fitness))) {
                    memcpy(particles[worst_i]->position,tp2,n_dimensions*sizeof(T));
                    particles[worst_i]->fitness = tp2_fitness;
                } else {
                    // Simplex komprimieren:
                    for (unsigned int i=0; i<n_dimensions; ++i) {
                        for (unsigned int j=0; j<n_points && j!=best_i; ++j) {
                            particles[j]->position[i] = particles[best_i]->position[i] 
                                                      + delta*(particles[j]->position[i]-particles[best_i]->position[i]);
                            // Keine min/max Pruefung beim Komprimieren (Grenzen koennen nicht ueberschritten werden)
                        }
                    }
                    for (unsigned int j=0; j<n_points && j!=best_i; ++j) {
                        particles[j]->fitness = fitness_function(particles[j]->position);
                    }
                }
            }
            // Indizes neu setzen:
            best_i = 0;
            worst_i = 0; second_worst_i = 0;
            for (unsigned int i=0; i<n_points; ++i) {
                if (search_min) {
                    if (particles[i]->fitness < particles[best_i]->fitness) {
                        best_i = i;
                    } else if (particles[i]->fitness > particles[worst_i]->fitness) {
                        second_worst_i = worst_i;
                        worst_i = i;
                    }
                } else {
                    if (particles[i]->fitness > particles[best_i]->fitness) {
                        best_i = i;
                    } else if (particles[i]->fitness < particles[worst_i]->fitness) {
                        second_worst_i = worst_i;
                        worst_i = i;
                    }
                }
            }
            if (worst_i == second_worst_i) {
                if (worst_i == 0) second_worst_i = 1;
                else second_worst_i = 0;
                for (unsigned int i=0; i<n_points && i!=worst_i; ++i) {
                    if (search_min) {
                        if (particles[i]->fitness > particles[second_worst_i]->fitness) second_worst_i = i;
                    } else {
                        if (particles[i]->fitness < particles[second_worst_i]->fitness) second_worst_i = i;
                    }
                }
            }
        }
    }
    delete center;
    delete tp;
    delete tp2;

    return total_best_particle->fitness;
}


template<class T,class U>
void NM_OPT<T,U>::generate_start_positions() {
    particles.push_back(new PARTICLE<T>(n_dimensions,total_best_particle->position));
    for (unsigned int i=1; i<n_points; ++i) {
        PARTICLE<T>* cp = new PARTICLE<T>(n_dimensions,total_best_particle->position);
        T tdim = 0.;
        T dim_dst = max_val[i]-min_val[i];
        do {
            tdim = min_val[i] + (rand()/T(RAND_MAX))*dim_dst;
        } while(fabs(cp->position[i-1]-tdim) < (dim_dst/20.));
        cp->position[i-1] = tdim;
        particles.push_back(cp);
    }
    if (!use_start_pos) { // Simplex verschieben, so dass start_pos im Mittelpunkt liegt
        PARTICLE<T>* cp = new PARTICLE<T>(n_dimensions);
        for (unsigned int i=0; i<n_dimensions; ++i) {
            cp->position[i] = 0.;
            for (unsigned int j=0; j<n_points; ++j) {
                cp->position[i] += particles[j]->position[i];
            }
            cp->position[i] /= n_points;
            cp->position[i] = total_best_particle->position[i] - cp->position[i];
        }
        for (unsigned int i=0; i<n_dimensions; ++i) {
            for (unsigned int j=0; j<n_points; ++j) {
                particles[j]->position[i] += cp->position[i];
            }
        }
        delete cp;
    }
}


//! ============================================================================


template<class T,class U>
CURVE_FIT<T,U>::CURVE_FIT(unsigned int const& n_dim,unsigned int const& n_facs,unsigned int const& n_data,
                          T const *const data,T const *const *const vars,U& mf,T const *const start_param,
                          T const *const min,T const *const max):n_dimensions(n_dim),n_factors(n_facs),
                          model_function(mf),n_data_points(n_data) {
    if (n_dim >= n_data) {
        cerr << "CURVE_FIT::CURVE_FIT -> warning: using only " << n_data << " data points to fit function with " << n_dimensions
             << " degrees of freedom.\n"
             << "                                 (you should use at least n_dimensions+1 data points)" << endl;
    }
    
    data_points = new T[n_data_points];
    memcpy(data_points,data,n_data_points*sizeof(T));

    data_vars = new T*[n_data_points];
    for (unsigned int i=0; i<n_data_points; ++i) {
        data_vars[i] = new T[n_dimensions];
        memcpy(data_vars[i],vars[i],n_dimensions*sizeof(T));
    }

    factors = new T[n_factors];
    if (start_param) {
        start_vals = new T[n_factors];
        memcpy(start_vals,start_param,n_factors*sizeof(T));
    } else start_vals = 0;

    min_val = new T[n_factors];
    if (min) memcpy(min_val,min,n_factors*sizeof(T));
    else for (unsigned int i=0; i<n_factors; ++i) min_val[i] = -999999999.;

    max_val = new T[n_factors];
    if (max) memcpy(max_val,max,n_factors*sizeof(T));
    else for (unsigned int i=0; i<n_factors; ++i) max_val[i] = 999999999.;

    bs_variances = new T[n_factors];
    bs_means = new T[n_factors];
    for (unsigned int i=0; i<n_factors; ++i) {
        bs_variances[i] = 0.;
        bs_means[i] = 0.;
    }
    bootstrap_calculated = false;
}


template<class T,class U>
CURVE_FIT<T,U>::~CURVE_FIT() {
    delete[] data_points;
    for (unsigned int i=0; i<n_data_points; ++i) delete[] data_vars[i];
    delete[] data_vars;
    delete[] factors;
    delete[] min_val;
    delete[] max_val;
    delete[] bs_variances;
    delete[] bs_means;
    if (start_vals) delete[] start_vals;
}


template<class T,class U>
T CURVE_FIT<T,U>::SQUARE_DIFF::operator()(T const *const vals) {
    T sdiff = 0.;
    for (unsigned int i=0; i<parent->n_data_points; ++i) {
        T sd = parent->data_points[i] - parent->model_function(parent->data_vars[i],vals);
        sd *= sd;
        sdiff += sd;
    }
    return sdiff;
}


template<class T,class U>
T CURVE_FIT<T,U>::fit(int const& quality) {
    used_quality = quality;    
    unsigned int n_particles = n_factors * 60;
    int n_global = 5;
    unsigned int max_iter = 1000;
    unsigned int no_imp = 30;
    if (quality == opt_gn_ns::low) {
        n_particles = n_factors * 25;
        n_global = 0;
        max_iter = 500;
        no_imp = 25;
    } else if (quality == opt_gn_ns::medium) {
        n_particles = n_factors * 30;
        n_global = 3;
        max_iter = 500;
        no_imp = 25;
    } else if (quality == opt_gn_ns::high) {
        n_particles = n_factors * 60;
        n_global = 5;
        max_iter = 1000;
        no_imp = 30;
    } else if (quality == opt_gn_ns::very_high) {
        n_particles = n_factors * 100;
        n_global = 10;
        max_iter = 5000;
        no_imp = 50;
    }

    SQUARE_DIFF sd(this);

    PS_OPT<T,SQUARE_DIFF> cps(n_factors,sd,start_vals,min_val,max_val,0,0,n_particles,max_iter,no_imp);

    cps.minimize(n_global);

    memcpy(factors,cps.get_best_position(),n_factors*sizeof(T));

    return sd(factors)/n_data_points;
}


template<class T,class U> T const *const CURVE_FIT<T,U>::get_factors() {return factors;}
template<class T,class U> T const CURVE_FIT<T,U>::get_error() {SQUARE_DIFF sd(this); return sd(factors)/(n_data_points - n_factors);}


template<class T,class U> 
void CURVE_FIT<T,U>::bootstrap(unsigned int const& n_simulations) {
    T** sim_factors = new T*[n_factors];
    for (unsigned int i=0; i<n_factors; ++i) {
        sim_factors[i] = new T[n_simulations];
    }
    T* residuals = new T[n_data_points];
    for (unsigned int j=0; j<n_data_points; ++j) {
        residuals[j] = data_points[j] - model_function(data_vars[j],factors);
    }
    for (unsigned int i=0; i<n_simulations; ++i) {
        // neues data_points samplen
        // und curve_fit machen --> variances updaten
        T* simdata = new T[n_data_points];
        
        /*
        // sample with replacement  (only with higher numbers of data points)        
        for (unsigned int j=0; j<n_data_points; ++j) {
            unsigned int z = rand() % n_data_points;
            simdata[j] = data_points[z];
        }
        */
        
        for (unsigned int j=0; j<n_data_points; ++j) {
            unsigned int z = rand() % n_data_points;
            simdata[j] = model_function(data_vars[j],factors) + residuals[z];
        }
        
        CURVE_FIT<T,U> simfit(n_dimensions,n_factors,n_data_points,simdata,data_vars,
                              model_function,start_vals,min_val,max_val);
        simfit.fit(used_quality);

        delete[] simdata;
        for (unsigned int j=0; j<n_factors; ++j) {
            sim_factors[j][i] = simfit.get_factors()[j]; 
        }
    }
    delete[] residuals;
    for (unsigned int i=0; i<n_factors; ++i) {
        bs_means[i] = 0.;
        for (unsigned int j=0; j<n_simulations; ++j) {
            bs_means[i] += sim_factors[i][j];
        }
        bs_means[i] /= n_simulations;
        bs_variances[i] = 0.;
        for (unsigned int j=0; j<n_simulations; ++j) {
            T sdiff = sim_factors[i][j] - bs_means[i];
            sdiff *= sdiff;
            bs_variances[i] += sdiff;
        }
        bs_variances[i] /= n_simulations - 1.;
    }
    for (unsigned int i=0; i<n_factors; ++i) delete[] sim_factors[i];
    delete[] sim_factors;
    
    bootstrap_calculated = true;
}


template<class T,class U> 
T const *const CURVE_FIT<T,U>::get_bootstrap_variances(unsigned int const& n_simulations) {
    if (bootstrap_calculated) return bs_variances;
    bootstrap(n_simulations);    
    return bs_variances;
}


template<class T,class U> 
T const *const CURVE_FIT<T,U>::get_bootstrap_means(unsigned int const& n_simulations) {
    if (bootstrap_calculated) return bs_means;
    bootstrap(n_simulations);    
    return bs_means;
}

#endif

