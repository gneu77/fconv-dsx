#ifndef __BITMATRIX_GN
#define __BITMATRIX_GN

#include<string.h>
#include<stdint.h>

using namespace std;

template<class T> //!interner Basistyp
class BITMATRIX {
    private:
        T* field;
        T field_size;
        T t_size;
        T teiler;
        T x_dim;
        T y_dim;
    public:
        BITMATRIX(uint32_t const& x,uint32_t const& y);
        ~BITMATRIX() {delete[] field;};
        
        inline void set_bit(T const& x,T const& y);
        inline T get_bit(T const& x,T const& y);
        
        inline void clear_bit(T const& x,T const& y);
        inline void flip_bit(T const& x,T const& y);
        inline void set_bit(int const& x,int const& y);
        inline T get_bit(int const& x,int const& y);
        
        inline void clear_bit(int const& x,int const& y);
        inline void flip_bit(int const& x,int const& y);
        inline void flip();
};

template<class T>
BITMATRIX<T>::BITMATRIX(uint32_t const& x,uint32_t const& y) {
    x_dim = x;
    y_dim = y;
    t_size = sizeof(T);
    teiler = t_size * 8;
    field_size = x * y;
    field_size /= teiler; field_size += 1;
    field = new T[field_size];
    memset(field,0,field_size*t_size);
}

template<class T>
void BITMATRIX<T>::set_bit(T const& x,T const& y) {
    T offset = y*x_dim + x;
    T index = offset / teiler;
    field[index] = field[index] | (0x01 << (offset % teiler));
}

template<class T>
T BITMATRIX<T>::get_bit(T const& x,T const& y) {
    T offset = y*x_dim + x;
    T index = offset / teiler;
    return (field[index] & (0x01 << (offset % teiler)));
}

template<class T>
void BITMATRIX<T>::clear_bit(T const& x,T const& y) {
    T offset = y*x_dim + x;
    T index = offset / teiler;
    field[index] = field[index] & (~(0x01 << (offset % teiler)));
}

template<class T>
void BITMATRIX<T>::flip_bit(T const& x,T const& y) {
    T offset = y*x_dim + x;
    T index = offset / teiler;
    field[index] = field[index] ^ (0x01 << (offset % teiler));
}

template<class T>
void BITMATRIX<T>::set_bit(int const& x,int const& y) {
    T offset = y*x_dim + x;
    T index = offset / teiler;
    field[index] = field[index] | (0x01 << (offset % teiler));
}

template<class T>
T BITMATRIX<T>::get_bit(int const& x,int const& y) {
    T offset = y*x_dim + x;
    T index = offset / teiler;
    return (field[index] & (0x01 << (offset % teiler)));
}

template<class T>
void BITMATRIX<T>::clear_bit(int const& x,int const& y) {
    T offset = y*x_dim + x;
    T index = offset / teiler;
    field[index] = field[index] & (~(0x01 << (offset % teiler)));
}

template<class T>
void BITMATRIX<T>::flip_bit(int const& x,int const& y) {
    T offset = y*x_dim + x;
    T index = offset / teiler;
    field[index] = field[index] ^ (0x01 << (offset % teiler));
}

template<class T>
void BITMATRIX<T>::flip() {
    for (T i=0; i<field_size; ++i) field[i] = ~(field[i]);
}

#endif

