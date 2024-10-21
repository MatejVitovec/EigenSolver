#ifndef FIELD_HPP
#define FIELD_HPP

#include <vector>
#include <memory>
#include <string>
#include <eigen3/Eigen/Dense>


template <int M, int N>
class Field
{
    typedef Eigen::Array<double, M, Eigen::Dynamic> FieldArray;

    public:
        Field() {}

        Field(size_t size) : data(FieldArray(M, size*N)) {}

        Field(FieldArray dataIn) : data(dataIn) {}

        const Eigen::Block<const FieldArray, M, N> operator[](int i) const
        {
            return Eigen::Block<const FieldArray, M, N>(data.derived(), 0, i*N);
        }

        Eigen::Block<FieldArray, M, N> operator[](int i)
        {
            return Eigen::Block<FieldArray, M, N>(data.derived(), 0, i*N);
        }

        FieldArray& getData()
        {
            return data;
        }

        void operator+=(const Field<M, N>& v);
        void operator-=(const Field<M, N>& v);
        void operator+=(const Eigen::MatrixBase<Field::FieldArray>& v);
        void operator-=(const Eigen::MatrixBase<Field::FieldArray>& v);

        template <int P, int Q> friend auto operator+(const Field<P, Q>& u, const Field<P, Q>& v);
        template <int P, int Q> friend auto operator-(const Field<P, Q>& u, const Field<P, Q>& v);
        template <int P, int Q> friend auto operator*(const Field<P, Q>& u, const Field<P, Q>& v);
        template <int P, int Q> friend auto operator/(const Field<P, Q>& u, const Field<P, Q>& v);
        template <int P, int Q> friend auto operator+(const Field<P, Q>& u, const Eigen::MatrixBase<Field::FieldArray>& v);
        template <int P, int Q> friend auto operator-(const Field<P, Q>& u, const Eigen::MatrixBase<Field::FieldArray>& v);
        template <int P, int Q> friend auto operator*(const Field<P, Q>& u, const Eigen::MatrixBase<Field::FieldArray>& v);
        template <int P, int Q> friend auto operator/(const Field<P, Q>& u, const Eigen::MatrixBase<Field::FieldArray>& v);
        template <int P, int Q> friend auto operator+(const Eigen::MatrixBase<Field::FieldArray>& u, const Field<P, Q>& v);
        template <int P, int Q> friend auto operator-(const Eigen::MatrixBase<Field::FieldArray>& u, const Field<P, Q>& v);
        template <int P, int Q> friend auto operator*(const Eigen::MatrixBase<Field::FieldArray>& u, const Field<P, Q>& v);
        template <int P, int Q> friend auto operator/(const Eigen::MatrixBase<Field::FieldArray>& u, const Field<P, Q>& v);

    protected:
        FieldArray data;
        
};

template <int M, int N>
void Field<M, N>::operator+=(const Field<M, N>& v)
{
    data += v.data;
}

template <int M, int N>
void Field<M, N>::operator-=(const Field<M, N>& v)
{
    data -= v.data;
}

template <int M, int N>
void Field<M, N>::operator+=(const Eigen::MatrixBase<Field::FieldArray>& v)
{
    data += v.data;
}

template <int M, int N>
void Field<M, N>::operator-=(const Eigen::MatrixBase<Field::FieldArray>& v)
{
    data -= v.data;
}

template <int M, int N>
auto operator+(const Field<M, N>& u, const Field<M, N>& v)
{
    return u.data + v.data;
}

template <int M, int N>
auto operator-(const Field<M, N>& u, const Field<M, N>& v)
{
    return u.data - v.data;
}

template <int M, int N>
auto operator*(const Field<M, N>& u, const Field<M, N>& v)
{
    return u.data * v.data;
}

template <int M, int N>
auto operator/(const Field<M, N>& u, const Field<M, N>& v)
{
    return u.data / v.data;
}

template <int M, int N>
auto operator+(const Field<M, N>& u, const Eigen::MatrixBase<Eigen::Array<double, M, Eigen::Dynamic>>& v)
{
    return u.data + v;
}

template <int M, int N>
auto operator-(const Field<M, N>& u, const Eigen::MatrixBase<Eigen::Array<double, M, Eigen::Dynamic>>& v)
{
    return u.data - v;
}

template <int M, int N>
auto operator*(const Field<M, N>& u, const Eigen::MatrixBase<Eigen::Array<double, M, Eigen::Dynamic>>& v)
{
    return u.data * v;
}

template <int M, int N>
auto operator/(const Field<M, N>& u, const Eigen::MatrixBase<Eigen::Array<double, M, Eigen::Dynamic>>& v)
{
    return u.data / v;
}

template <int M, int N>
auto operator+(const Eigen::MatrixBase<Eigen::Array<double, M, Eigen::Dynamic>>& u, const Field<M, N>& v)
{
    return u + v.data;
}

template <int M, int N>
auto operator-(const Eigen::MatrixBase<Eigen::Array<double, M, Eigen::Dynamic>>& u, const Field<M, N>& v)
{
    return u - v.data;
}

template <int M, int N>
auto operator*(const Eigen::MatrixBase<Eigen::Array<double, M, Eigen::Dynamic>>& u, const Field<M, N>& v)
{
    return u * v.data;
}

template <int M, int N>
auto operator/(const Eigen::MatrixBase<Eigen::Array<double, M, Eigen::Dynamic>>& u, const Field<M, N>& v)
{
    return u / v.data;
}




#endif // FIELD_HPP