#ifndef MATH_H
#define MATH_H

#include <cstdio>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <vector>

namespace ZMATH
{

#define PI                      3.14159265358979323
#define DEG2RAD(theta)          (0.01745329251994329 * (theta))
#define RAD2DEG(theta)          (57.2957795130823208 * (theta))
#define WGS_84_RADIUS_EQUATOR   6378137.0
#define WGS_84_RADIUS_POLAR     6356752.3142

	typedef unsigned char           byte;
	typedef unsigned short          ushort;
	typedef unsigned int            uint;
	typedef unsigned long           ulong;

	template<class T>
	inline T tmin(T a, T b)
	{
		return a < b ? a : b;
	}

	template<class T>
	inline T tmax(T a, T b)
	{
		return a > b ? a : b;
	}

	//-------------------- 2维坐标类
	template <typename T>
	class tvector2
	{
	public:
		typedef T value_type;
		typedef std::size_t size_type;
		typedef tvector2<T> type;

		//-------------------- point
		value_type x;
		value_type y;
		size_type length() const
		{
			return 2;
		}

		//-------------------- 构造函数
		tvector2() :
			x(value_type(0)),
			y(value_type(0)) { }
		tvector2(type const & v) :
			x(v.x),
			y(v.y) { }
		tvector2(value_type const & s) :
			x(s),
			y(s) { }

		tvector2(value_type const & s1, value_type const & s2) :
			x(s1),
			y(s2) { }

		template <typename U>
		tvector2(U const & x) :
			x(value_type(x)),
			y(value_type(x)) { }

		template <typename U, typename V>
		tvector2(U const & a, V const & b) :
			x(value_type(a)),
			y(value_type(b)) { }

		template <typename U>
		tvector2(tvector2<U> const & v) :
			x(value_type(v.x)),
			y(value_type(v.y)) { }

		//-------------------- 运算符重载
		value_type & operator[](size_type i)
		{
			assert(i < this->length());
			return (&x)[i];
		}

		value_type const & operator[](size_type i) const
		{
			assert(i < this->length());
			return (&x)[i];
		}

		type & operator = (type const & v)
		{
			this->x = v.x;
			this->y = v.y;
			return *this;
		}

		template <typename U>
		type & operator = (tvector2<U> const & v)
		{
			this->x = T(v.x);
			this->y = T(v.y);
			return *this;
		}

		template <typename U>
		type & operator += (U const & s)
		{
			this->x += T(s);
			this->y += T(s);
			return *this;
		}

		template <typename U>
		type & operator += (tvector2<U> const & v)
		{
			this->x += T(v.x);
			this->y += T(v.y);
			return *this;
		}


		template <typename U>
		type & operator -= (U const & s)
		{
			this->x -= T(s);
			this->y -= T(s);
			return *this;
		}

		template <typename U>
		type & operator -= (tvector2<U> const & v)
		{
			this->x -= T(v.x);
			this->y -= T(v.y);
			return *this;
		}

		template <typename U>
		type & operator *= (U s)
		{
			this->x *= T(s);
			this->y *= T(s);
			return *this;
		}

		template <typename U>
		type & operator *= (tvector2<U> const & v)
		{
			this->x *= T(v.x);
			this->y *= T(v.y);
			return *this;
		}

		template <typename U>
		type & operator /= (U s)
		{
			this->x /= T(s);
			this->y /= T(s);
			return *this;
		}

		template <typename U>
		type & operator /= (tvector2<U> const & v)
		{
			this->x /= T(v.x);
			this->y /= T(v.y);
			return *this;
		}

		type & operator++()
		{
			++ this->x;
			++ this->y;
			return *this;
		}
		type & operator--()
		{
			-- this->x;
			-- this->y;
			return *this;
		}

		friend bool operator == (type const & v1, type const & v2)
		{
			return (v1.x == v2.x) && (v1.y == v2.y);
		}

		friend bool operator != (type const & v1, type const & v2)
		{
			return (v1.x != v2.x) || (v1.y != v2.y);
		}

		type operator + (T const & s)
		{
			return type(
				x + T(s),
				y + T(s));
		}

		friend type operator + (T const & s, type const & v)
		{
			return v + s;
		}

		type operator + (type const & v)
		{
			return type(
				x + T(v.x),
				y + T(v.y));
		}

		type operator - (T const & s)
		{
			return type(
				x - T(s),
				y - T(s));
		}

		friend type operator - (T const & s, type const & v)
		{
			return type(
				T(s) - v.x,
				T(s) - v.y);
		}

		type operator - (type const & v)
		{
			return type(
				x - T(v.x),
				y - T(v.y));
		}

		type operator * (T const & s)
		{
			return type(
				x * T(s),
				y * T(s));
		}

		friend type operator * (T const & s, type const & v)
		{
			return v * s;
		}

		type operator * (type const & v)
		{
			return type(
				x * T(v.x),
				y * T(v.y));
		}

		type operator / (T const & s)
		{
			return type(
				x / T(s),
				y / T(s));
		}

		friend type operator / (T const & s, type const & v)
		{
			return type(
				T(s) / v.x,
				T(s) / v.y);
		}

		type operator / (type const & v)
		{
			return  type(
				x / T(v.x),
				y / T(v.y)
			);
		}

		type operator - ()
		{
			return  type(
				-x,
				-y
			);
		}

		type operator ++ (int)
		{
			return  *this + T(1);
		}

		type operator -- (int)
		{
			return *this - T(1);
		}

	};

	template<typename T>
	tvector2<T> rotate(tvector2<T> const & v, T angle)
	{
		tvector2<T> res;
		T const c(cos(DEG2RAD(angle)));
		T const s(sin(DEG2RAD(angle)));
		res.x = v.x * c - v.y * s;
		res.y = v.x * s + v.y * c;
		return res;
	}

	//-------------------- 3维坐标类
	template<typename T>
	class tvector3
	{
	public:
		typedef T value_type;
		typedef std::size_t size_type;
		typedef tvector3<T> type;

		//-------------------- point
		value_type x;
		value_type y;
		value_type z;
		size_type length() const
		{
			return 3;
		}

		//-------------------- 构造函数
		inline tvector3() :
			x(value_type(0)),
			y(value_type(0)),
			z(value_type(0)) { }

		inline tvector3(type const & v) :
			x(v.x),
			y(v.y),
			z(v.z) { }

		inline tvector3(value_type s) :
			x(s),
			y(s),
			z(s) { }

		inline tvector3(value_type s0, value_type s1, value_type s2) :
			x(s0),
			y(s1),
			z(s2) { }

		template <typename U>
		tvector3(U s) :
			x(value_type(s)),
			y(value_type(s)),
			z(value_type(s)) { }

		template <typename A, typename B, typename C>
		tvector3(A x, B y, C z) :
			x(value_type(x)),
			y(value_type(y)),
			z(value_type(z)) { }

		template <typename A, typename B>
		tvector3(tvector3<A> const& v, B s) :
			x(value_type(v.x)),
			y(value_type(v.y)),
			z(value_type(s)) { }
		template <typename A, typename B>
		tvector3(A s, tvector3<B> const& v) :
			x(value_type(s)),
			y(value_type(v.x)),
			z(value_type(v.y)) { }

		template <typename U>
		tvector3(tvector3<U> const & v) :
			x(value_type(v.x)),
			y(value_type(v.y)),
			z(value_type(v.z)) { }

		//-------------------- 运算符重载
		value_type & operator[](size_type i)
		{
			assert(i < this->length());
			return (&x)[i];
		}

		value_type const & operator[](size_type i) const
		{
			assert(i < this->length());
			return (&x)[i];
		}

		type & operator = (type const & v)
		{
			this->x = v.x;
			this->y = v.y;
			this->z = v.z;
			return *this;
		}
		template <typename U>
		type & operator = (type const & v)
		{
			this->x = T(v.x);
			this->y = T(v.y);
			this->z = T(v.z);
			return *this;
		}

		template <typename U>
		type & operator += (U const & s)
		{
			this->x += T(s);
			this->y += T(s);
			this->z += T(s);
			return *this;
		}

		template <typename U>
		type & operator += (type const & v)
		{
			this->x += T(v.x);
			this->y += T(v.y);
			this->z += T(v.z);
			return *this;
		}

		template <typename U>
		type & operator -= (U const & s)
		{
			this->x -= T(s);
			this->y -= T(s);
			this->z -= T(s);
			return *this;
		}

		template <typename U>
		type & operator -= (type const & v)
		{
			this->x -= T(v.x);
			this->y -= T(v.y);
			this->z -= T(v.z);
			return *this;
		}
		template <typename U>
		type & operator *= (U const & s)
		{
			this->x *= T(s);
			this->y *= T(s);
			this->z *= T(s);
			return *this;
		}

		template <typename U>
		type & operator *= (type const & v)
		{
			this->x *= T(v.x);
			this->y *= T(v.y);
			this->z *= T(v.z);
			return *this;
		}

		template <typename U>
		type & operator /= (U const & s)
		{
			this->x /= T(s);
			this->y /= T(s);
			this->z /= T(s);
			return *this;
		}

		template <typename U>
		type & operator /= (type const & v)
		{
			this->x /= T(v.x);
			this->y /= T(v.y);
			this->z /= T(v.z);
			return *this;
		}

		type & operator ++ ()
		{
			++this->x;
			++this->y;
			++this->z;
			return *this;
		}
		type & operator -- ()
		{
			--this->x;
			--this->y;
			--this->z;
			return *this;
		}

		void makeFloor(const type& cmp)
		{
			if (cmp.x < x) x = cmp.x;
			if (cmp.y < y) y = cmp.y;
			if (cmp.z < z) z = cmp.z;
		}
		void makeCeil(const type& cmp)
		{
			if (cmp.x > x) x = cmp.x;
			if (cmp.y > y) y = cmp.y;
			if (cmp.z > z) z = cmp.z;
		}
		T lengthf() const
		{
			return (T)sqrtf(x * x + y * y + z * z);
		}

		bool operator > (const type & r) const
		{
			return x > r.x && y > r.y && z > r.z;
		}

		bool operator < (const type & r) const
		{
			return x < r.x && y < r.y && z < r.z;
		}

		bool operator == (type const & v2) const
		{
			return (x == v2.x) && (y == v2.y) && (z == v2.z);
		}

		bool operator != (type const & v2) const
		{
			return (x != v2.x) || (y != v2.y) || (z != v2.z);
		}

		type operator + (T const & s) const
		{
			return type(
				x + T(s),
				y + T(s),
				z + T(s));
		}

		friend type operator + (T const & s, type const & v)
		{
			return v + s;
		}

		type operator + (type const & v2) const
		{
			return type(
				x + T(v2.x),
				y + T(v2.y),
				z + T(v2.z));
		}

		type operator - (T const & s) const
		{
			return type(
				x - T(s),
				y - T(s),
				z - T(s));
		}

		friend type operator - (T const & s, type const & v)
		{
			return type(
				T(s) - v.x,
				T(s) - v.y,
				T(s) - v.z);
		}

		type operator - (type const & v2) const
		{
			return type(
				x - T(v2.x),
				y - T(v2.y),
				z - T(v2.z));
		}

		type operator * (T const & s) const
		{
			return type(
				x * T(s),
				y * T(s),
				z * T(s));
		}

		friend type operator * (T const & s, type const & v)
		{
			return v * s;
		}

		type operator * (type const & v2) const
		{
			return type(
				x * T(v2.x),
				y * T(v2.y),
				z * T(v2.z));
		}

		type operator / (T const & s) const
		{
			return type(
				x / T(s),
				y / T(s),
				z / T(s));
		}

		friend type operator / (T const & s, type const & v)
		{
			return type(
				T(s) / v.x,
				T(s) / v.y,
				T(s) / v.z);
		}

		type operator / (type const & v2) const
		{
			return type(
				x / T(v2.x),
				y / T(v2.y),
				z / T(v2.z));
		}

		type operator - () const
		{
			return type(-x, -y, -z);
		}

		type operator ++ (int) const
		{
			return type(
				x + T(1),
				y + T(1),
				z + T(1));
		}

		type operator -- (int) const
		{
			return type(
				x - T(1),
				y - T(1),
				z - T(1));
		}

		//3D空间旋转指定的角度，沿着垂直于X轴的方向顺指针旋转
		friend type rotateX(const type & v, T angle)
		{
			type res(v);
			T c = cos(T(DEG2RAD(angle)));
			T s = sin(T(DEG2RAD(angle)));

			res.y = v.y * c - v.z * s;
			res.z = v.y * s + v.z * c;
			return res;
		}

		friend type rotateY(type const & v, T angle)
		{
			type res = v;
			T c = cos(T(DEG2RAD(angle)));
			T s = sin(T(DEG2RAD(angle)));

			res.x = v.x * c + v.z * s;
			res.z = -v.x * s + v.z * c;
			return res;
		}

		friend type rotateZ(type const & v, T angle)
		{
			type res = v;
			T c = cos(DEG2RAD(angle));
			T s = sin(DEG2RAD(angle));

			res.x = v.x * c - v.y * s;
			res.y = v.x * s + v.y * c;
			return res;
		}
	};

	template <typename T>
	tvector3<T> rotateX(const tvector3<T>& v, T angle)
	{
		tvector3<T> res(v);
		T c = cos(T(DEG2RAD(angle)));
		T s = sin(T(DEG2RAD(angle)));

		res.y = v.y * c - v.z * s;
		res.z = v.y * s + v.z * c;
		return res;
	}

	template <typename T>
	tvector3<T> rotateY(tvector3<T> const & v, T angle)
	{
		tvector3<T> res = v;

		T c = cos(T(DEG2RAD(angle)));
		T s = sin(T(DEG2RAD(angle)));

		res.x = v.x * c + v.z * s;
		res.z = -v.x * s + v.z * c;
		return res;
	}

	template <typename T>
	tvector3<T> rotateZ(tvector3<T> const & v, T angle)
	{

		tvector3<T> res = v;

		T c = cos(DEG2RAD(angle));

		T s = sin(DEG2RAD(angle));

		res.x = v.x * c - v.y * s;
		res.y = v.x * s + v.y * c;
		return res;
	}

	//向量夹角
	/**
	*   两个向量的夹角
	*   A·B = |A|*|B|*cos(θ)
	*   cos(θ) = A·B / |A|*|B|
	*/
	template<typename T>
	inline bool _isnan(T t)
	{
		return t == t;
	}

	template <typename T>
	inline T acosEx(T val)
	{
		if (T(-1.0f) < val)
		{
			if (val < 1.0f)
				return T(acos(val));
			else
				return T(0);
		}
		else
		{
			return T(PI);
		}
	}

	template<typename T>
	T angleBetweenVector(const tvector3<T> & a, const tvector3<T> & b)
	{
#define Mag(V) (sqrtf(V.x*V.x + V.y*V.y + V.z*V.z))
		T dotProduct = dot(a, b);
		T vectorsMagnitude = Mag(a) * Mag(b);
		T angle = acos(dotProduct / vectorsMagnitude);
		T result = angle * T(RAD2DEG);
		if (_isnan(result))
		{
			return T(0);
		}
		else
		{
			return result;
		}
	}

	template <typename T>
	T angleBetweenVector(const tvector2<T> & a, const tvector2<T> & b)
	{
#define Mag2D(V) (sqrtf(V.x*V.x + V.y*V.y))
		T dotProduct = dot(a, b);
		T vectorsMagnitude = Mag2D(a) * Mag2D(b);
		T angle = acos(dotProduct / vectorsMagnitude);
		T result = angle * T(RAD2DEG);
		if (_isnan(result))
		{
			return T(0);
		}
		else
		{
			return result;
		}
	}

	////判断点与多边形位置
	///**
	//*	点在多边形里
	//*   当且仅当点与边的夹角之和 == 360
	//*/
	//template<typename T>
	//bool insidePolygon(const tvector3<T>& point, const tvector3<T> polygon[], size_t count)
	//{
	//	tvector3<T> vA, vB;
	//	T angle = T(0.0);
	//	for (size_t i = 0; i < count; ++i)
	//	{
	//		vA = polygon[i] - point;
	//		vB = polygon[(i + 1) % count] - point;
	//		angle += angleBetweenVector(vA, vB);
	//	}
	//	if (abs(angle - 360) >= 0.5f)
	//	{
	//		return true;
	//	}
	//	return false;
	//}

	//template<typename T>
	//bool insidePolygon(const tvector2<T>& point, const tvector2<T> polygon[], size_t count)
	//{
	//	T angle = T(0.0);
	//	tvector2<T> vA, vB;
	//	for (size_t i = 0; i < count; ++i)
	//	{
	//		vA = polygon[i] - point;
	//		vB = polygon[(i + 1) % count] - point;
	//		tvector3<T> a(vA.x, vA.y, 0);
	//		tvector3<T> b(vB.x, vB.y, 0);
	//		angle += angleBetweenVector<T>(a, b);
	//	}
	//	if (abs(angle - 360) >= 0.5f)
	//	{
	//		return true;
	//	}
	//	return false;
	//}

	//template<typename T>
	//bool pointinTriangle(tvector3<T> A, tvector3<T> B, tvector3<T> C, tvector3<T> P)
	//{
	//	tvector3<T> v0 = C - A;
	//	tvector3<T> v1 = B - A;
	//	tvector3<T> v2 = P - A;

	//	float dot00 = dot(v0, v0);
	//	float dot01 = dot(v0, v1);
	//	float dot02 = dot(v0, v2);
	//	float dot11 = dot(v1, v1);
	//	float dot12 = dot(v1, v2);

	//	float inverDeno = 1 / (dot00 * dot11 - dot01 * dot01);

	//	float u = (dot11 * dot02 - dot01 * dot12) * inverDeno;
	//	if (u < 0 || u > 1) // if u out of range, return directly
	//	{
	//		return false;
	//	}

	//	float v = (dot00 * dot12 - dot01 * dot02) * inverDeno;
	//	if (v < 0 || v > 1) // if v out of range, return directly
	//	{
	//		return false;
	//	}

	//	return u + v <= 1;
	//}

	//template<typename T>
	//bool pointinTriangle(tvector2<T> A, tvector2<T> B, tvector2<T> C, tvector2<T> P)
	//{
	//	return pointinTriangle(
	//		tvector3<T>(A.x, A.y, 0),
	//		tvector3<T>(B.x, B.y, 0),
	//		tvector3<T>(C.x, C.y, 0),
	//		tvector3<T>(P.x, P.y, 0));
	//}

	//-------------------- 4维坐标类
	template <typename T>
	class tvector4
	{
	public:
		typedef T value_type;
		typedef std::size_t size_type;
		typedef tvector4<T> type;

		//-------------------- point
		value_type x;
		value_type y;
		value_type z;
		value_type w;

		size_type length() const
		{
			return 4;
		}

		//-------------------- 构造函数
		tvector4() :
			x(value_type(0)),
			y(value_type(0)),
			z(value_type(0)),
			w(value_type(0)) { }
		tvector4(tvector3<T> const& v, T s) :
			x(v.x),
			y(v.y),
			z(v.z),
			w(s) { }

		tvector4(T s) :
			x(s),
			y(s),
			z(s),
			w(s) { }

		tvector4(type const & v) :
			x(v.x),
			y(v.y),
			z(v.z),
			w(v.w) { }

		template <typename A, typename B>
		tvector4(tvector3<A> const & v, B s) :
			x(value_type(v.x)),
			y(value_type(v.y)),
			z(value_type(v.z)),
			w(value_type(s)) { }

		template <typename A, typename B>
		tvector4(A s, tvector3<B> const & v) :
			x(value_type(s)),
			y(value_type(v.x)),
			z(value_type(v.y)),
			w(value_type(v.z)) { }

		template<typename U>
		tvector4(tvector4<U> const & v) :
			x(value_type(v.x)),
			y(value_type(v.y)),
			z(value_type(v.z)),
			w(value_type(v.w)) { }

		tvector4
		(
			value_type s1,
			value_type s2,
			value_type s3,
			value_type s4
		) :
			x(s1),
			y(s2),
			z(s3),
			w(s4) { }

		//-------------------- 运算符重载
		value_type & operator[](size_type i)
		{
			assert(i < this->length());
			return (&x)[i];
		}

		value_type const & operator[](size_type i) const
		{
			assert(i < this->length());
			return (&x)[i];
		}

		type & operator=(type const & v)
		{
			this->x = v.x;
			this->y = v.y;
			this->z = v.z;
			this->w = v.w;
			return *this;
		}

		template <typename U>
		type & operator= (tvector4<U> const & v)
		{
			this->x = T(v.x);
			this->y = T(v.y);
			this->z = T(v.z);
			this->w = T(v.w);
			return *this;
		}

		template <typename U>
		type & operator+=(U const & s)
		{
			this->x += T(s);
			this->y += T(s);
			this->z += T(s);
			this->w += T(s);
			return *this;
		}
		template <typename U>
		type & operator+=(tvector4<U> const & v)
		{
			this->x += T(v.x);
			this->y += T(v.y);
			this->z += T(v.z);
			this->w += T(v.w);
			return *this;
		}

		template <typename U>
		type & operator-=(U const & s)
		{
			this->x -= T(s);
			this->y -= T(s);
			this->z -= T(s);
			this->w -= T(s);
			return *this;
		}

		template <typename U>
		type & operator-=(tvector4<U> const & v)
		{
			this->x -= T(v.x);
			this->y -= T(v.y);
			this->z -= T(v.z);
			this->w -= T(v.w);
			return *this;
		}
		template <typename U>
		type & operator*= (U const & s)
		{
			this->x *= T(s);
			this->y *= T(s);
			this->z *= T(s);
			this->w *= T(s);
			return *this;
		}

		template <typename U>
		type & operator*=(tvector4<U> const & v)
		{
			this->x *= T(v.x);
			this->y *= T(v.y);
			this->z *= T(v.z);
			this->w *= T(v.w);
			return *this;
		}
		template <typename U>
		type & operator/=(U const & s)
		{
			this->x /= T(s);
			this->y /= T(s);
			this->z /= T(s);
			this->w /= T(s);
			return *this;
		}

		template <typename U>
		type & operator/=(tvector4<U> const & v)
		{
			this->x /= T(v.x);
			this->y /= T(v.y);
			this->z /= T(v.z);
			this->w /= T(v.w);
			return *this;
		}

		type & operator++()
		{
			++this->x;
			++this->y;
			++this->z;
			++this->w;
			return *this;
		}

		type & operator--()
		{
			--this->x;
			--this->y;
			--this->z;
			--this->w;
			return *this;
		}

		type operator + (T const & s) const
		{
			return type(
				x + s,
				y + s,
				z + s,
				w + s);
		}

		friend type operator + (T const & s, type const & v)
		{
			return v + s;
		}

		type operator + (type const & v2) const
		{
			return type(
				x + v2.x,
				y + v2.y,
				z + v2.z,
				w + v2.w);
		}

		type operator - (T const & s) const
		{
			return type(
				x - s,
				y - s,
				z - s,
				w - s);
		}

		friend type operator - (T const & s, type const & v)
		{
			return type(
				s - v.x,
				s - v.y,
				s - v.z,
				s - v.w);
		}

		type operator - (type const & v2) const
		{
			return type(
				x - v2.x,
				y - v2.y,
				z - v2.z,
				w - v2.w);
		}

		type operator * (T const & s) const
		{
			return type(
				x * s,
				y * s,
				z * s,
				w * s);
		}

		friend type operator * (T const & s, type const & v)
		{
			return v * s;
		}

		type operator * (type const & v2)
		{
			return type(
				x * v2.x,
				y * v2.y,
				z * v2.z,
				w * v2.w);
		}

		type operator / (T const & s) const
		{
			return type(
				x / s,
				y / s,
				z / s,
				w / s);
		}

		friend type operator / (T const & s, type const & v)
		{
			return type(
				s / v.x,
				s / v.y,
				s / v.z,
				s / v.w);
		}

		type operator / (type const & v2) const
		{
			return type(
				x / v2.x,
				y / v2.y,
				z / v2.z,
				w / v2.w);
		}

		type operator - () const
		{
			return type(
				-x,
				-y,
				-z,
				-w);
		}

		bool operator == (type const & v2) const
		{
			return (x == v2.x) && (y == v2.y) && (z == v2.z) && (w == v2.w);
		}

		bool operator != (type const & v2) const
		{
			return (x != v2.x) || (y != v2.y) || (z != v2.z) || (w != v2.w);
		}

		//3D空间旋转指定的角度，沿着垂直于X轴的方向顺指针旋转
		friend type rotateX(const type & v, T angle)
		{
			type res(v);
			T c = cos(T(DEG2RAD(angle)));
			T s = sin(T(DEG2RAD(angle)));

			res.y = v.y * c - v.z * s;
			res.z = v.y * s + v.z * c;
			return res;
		}

		friend type rotateY(type const & v, T angle)
		{
			type res = v;
			T c = cos(T(DEG2RAD(angle)));
			T s = sin(T(DEG2RAD(angle)));

			res.x = v.x * c + v.z * s;
			res.z = -v.x * s + v.z * c;
			return res;
		}

		friend type rotateZ(type const & v, T angle)
		{
			type res = v;
			T c = cos(DEG2RAD(angle));
			T s = sin(DEG2RAD(angle));

			res.x = v.x * c - v.y * s;
			res.y = v.x * s + v.y * c;
			return res;
		}
	};

	//-------------------- 2x2 Matrix类（列优先）
	template <typename T>
	class tmatrix2x2
	{
	public:
		typedef T				value_type;
		typedef std::size_t		size_type;
		typedef tvector2<T>	col_type;
		typedef tvector2<T>	row_type;
		typedef tmatrix2x2<T>	type;
		typedef tmatrix2x2<T>	transpose_type;
	private:
		col_type value[2];

	public:
		size_type length() const
		{
			return 2;
		}
		static size_type col_size()
		{
			return 2;
		}

		static size_type row_size()
		{
			return 2;
		}

		col_type & operator[](size_type i)
		{
			assert(i < this->length());
			return this->value[i];
		}
		col_type const & operator[](size_type i) const
		{
			assert(i < this->length());
			return this->value[i];
		}

		//-------------------- 求逆
		tmatrix2x2<T> _inverse() const
		{
			value_type Determinant = this->value[0][0] * this->value[1][1] - this->value[1][0] * this->value[0][1];

			tmat2x2<T> Inverse(
				+this->value[1][1] / Determinant,
				-this->value[0][1] / Determinant,
				-this->value[1][0] / Determinant,
				+this->value[0][0] / Determinant);
			return Inverse;
		}

		//-------------------- 构造函数
		tmatrix2x2()
		{
			this->value[0] = col_type(1, 0);
			this->value[1] = col_type(0, 1);
		}

		tmatrix2x2(type const & m)
		{
			this->value[0] = m.value[0];
			this->value[1] = m.value[1];
		}

		tmatrix2x2(value_type s)
		{
			this->value[0] = col_type(s, T(0));
			this->value[1] = col_type(T(0), s);
		}

		tmatrix2x2(value_type x0, value_type y0, value_type x1, value_type y1)
		{
			this->value[0] = col_type(x0, y0);
			this->value[1] = col_type(x1, y1);
		}

		tmatrix2x2(col_type const & v0, col_type const & v1)
		{
			this->value[0] = v0;
			this->value[1] = v1;
		}

		template <typename U>
		tmatrix2x2(U s)
		{
			this->value[0] = tvector2<T>(value_type(s), T(0));
			this->value[1] = tvector2<T>(T(0), value_type(s));
		}

		template <typename X1, typename Y1, typename X2, typename Y2>
		tmatrix2x2(X1 x1, Y1 y1, X2 x2, Y2 y2)
		{
			this->value[0] = col_type(value_type(x1), value_type(y1));
			this->value[1] = col_type(value_type(x2), value_type(y2));
		}

		template <typename V1, typename V2>
		tmatrix2x2
		(
			tvector2<V1> const & v1,
			tvector2<V2> const & v2
		)
		{
			this->value[0] = col_type(v1);
			this->value[1] = col_type(v2);
		}

		template <typename U>
		tmatrix2x2(tmatrix2x2<U> const & m)
		{
			this->value[0] = col_type(m[0]);
			this->value[1] = col_type(m[1]);
		}

		//-------------------- 运算符重载
		type & operator = (type const & m)
		{
			this->value[0] = m[0];
			this->value[1] = m[1];
			return *this;
		}

		template <typename U>
		type & operator = (tmatrix2x2<U> const & m)
		{
			this->value[0] = m[0];
			this->value[1] = m[1];
			return *this;
		}

		template <typename U>
		type & operator += (U const & s)
		{
			this->value[0] += s;
			this->value[1] += s;
			return *this;
		}

		template <typename U>
		type & operator += (tmatrix2x2<U> const & m)
		{
			this->value[0] += m[0];
			this->value[1] += m[1];
			return *this;
		}

		template <typename U>
		type & operator -= (U const & s)
		{
			this->value[0] -= s;
			this->value[1] -= s;
			return *this;
		}

		template <typename U>
		type & operator -= (tmatrix2x2<U> const & m)
		{
			this->value[0] -= m[0];
			this->value[1] -= m[1];
			return *this;
		}

		template <typename U>
		type & operator *= (U const & s)
		{
			this->value[0] *= s;
			this->value[1] *= s;
			return *this;
		}

		template <typename U>
		type & operator *= (tmatrix2x2<U> const & m)
		{
			return (*this = *this * m);
		}

		template <typename U>
		type & operator /= (U const & s)
		{
			this->value[0] /= s;
			this->value[1] /= s;
			return *this;
		}

		template <typename U>
		type & operator /= (tmatrix2x2<U> const & m)
		{
			return (*this = *this / m);
		}

		type & operator ++ ()
		{
			++this->value[0];
			++this->value[1];
			return *this;
		}

		type & operator -- ()
		{
			--this->value[0];
			--this->value[1];
			return *this;
		};

		friend tmatrix2x2<T> rotate(T angle)
		{
			T c = cos(DEG2RAD(angle));
			T s = sin(DEG2RAD(angle));
			return tmatrix2x2<T>(c, s, -s, c);
		}

		friend tmatrix2x2<T> operator + (tmatrix2x2<T> const & m, T const & s)
		{
			return tmatrix2x2<T>(m[0] + s, m[1] + s);
		}

		friend tmatrix2x2<T> operator + (T const & s, tmatrix2x2<T> const & m)
		{
			return tmatrix2x2<T>(m[0] + s, m[1] + s);
		}

		friend tmatrix2x2<T> operator + (tmatrix2x2<T> const & m1, tmatrix2x2<T> const & m2)
		{
			return tmatrix2x2<T>(m1[0] + m2[0], m1[1] + m2[1]);
		}

		friend tmatrix2x2<T> operator - (tmatrix2x2<T> const & m, T const & s)
		{
			return tmatrix2x2<T>(m[0] - s, m[1] - s);
		}

		friend tmatrix2x2<T> operator - (T const & s, tmatrix2x2<T> const & m)
		{
			return tmatrix2x2<T>(s - m[0], s - m[1]);
		}

		friend tmatrix2x2<T> operator - (tmatrix2x2<T> const & m1, tmatrix2x2<T> const & m2)
		{
			return tmatrix2x2<T>(m1[0] - m2[0], m1[1] - m2[1]);
		}

		friend tmatrix2x2<T> operator * (tmatrix2x2<T> const & m, T  const & s)
		{
			return tmatrix2x2<T>(m[0] * s, m[1] * s);
		}

		friend tmatrix2x2<T> operator * (T const & s, tmatrix2x2<T> const & m)
		{
			return tmatrix2x2<T>(m[0] * s, m[1] * s);
		}

		friend tvector2<T> operator * (tmatrix2x2<T> const & m, tvector2<T> const & v)
		{
			return tvector2<T>(
				m[0][0] * v.x + m[1][0] * v.y,
				m[0][1] * v.x + m[1][1] * v.y);
		}

		friend tvector2<T> operator * (tvector2<T> const & v, tmatrix2x2<T> const & m)
		{
			return  tvector2<T>(
				v.x * m[0][0] + v.y * m[0][1],
				v.x * m[1][0] + v.y * m[1][1]);
		}

		friend tmatrix2x2<T> operator * (tmatrix2x2<T> const & m1, tmatrix2x2<T> const & m2)
		{
			return tmatrix2x2<T>(
				m1[0][0] * m2[0][0] + m1[1][0] * m2[0][1],
				m1[0][1] * m2[0][0] + m1[1][1] * m2[0][1],
				m1[0][0] * m2[1][0] + m1[1][0] * m2[1][1],
				m1[0][1] * m2[1][0] + m1[1][1] * m2[1][1]);
		}
	};

	//-------------------- 3x3 Matrix类（列优先）
	template <typename T>
	class tmatrix3x3
	{
	public:
		typedef T				value_type;
		typedef std::size_t		size_type;
		typedef tvector3<T>	col_type;
		typedef tvector3<T>	row_type;
		typedef tmatrix3x3<T>	type;
		typedef tmatrix3x3<T>	transpose_type;

	private:
		col_type value[3];

	public:
		size_type length() const
		{
			return 3;
		}
		size_type col_size()
		{
			return 3;
		}
		size_type row_size()
		{
			return 3;
		}

		//-------------------- 构造函数
		tmatrix3x3()
		{
			this->value[0] = col_type(T(1), T(0), T(0));
			this->value[1] = col_type(T(0), T(1), T(0));
			this->value[2] = col_type(T(0), T(0), T(1));
		}

		tmatrix3x3(type const & m)
		{
			this->value[0] = m.value[0];
			this->value[1] = m.value[1];
			this->value[2] = m.value[2];
		}

		tmatrix3x3(value_type const & s)
		{
			this->value[0] = col_type(s, T(0), T(0));
			this->value[1] = col_type(T(0), s, T(0));
			this->value[2] = col_type(T(0), T(0), s);
		}

		tmatrix3x3
		(
			value_type const & x0, value_type const & y0, value_type const & z0,
			value_type const & x1, value_type const & y1, value_type const & z1,
			value_type const & x2, value_type const & y2, value_type const & z2
		)
		{
			this->value[0] = col_type(x0, y0, z0);
			this->value[1] = col_type(x1, y1, z1);
			this->value[2] = col_type(x2, y2, z2);
		}

		tmatrix3x3
		(
			col_type const & v0,
			col_type const & v1,
			col_type const & v2
		)
		{
			this->value[0] = v0;
			this->value[1] = v1;
			this->value[2] = v2;
		}

		template <typename U>
		tmatrix3x3(U const & s)
		{
			this->value[0] = tvector3<T>(T(s), T(0), T(0));
			this->value[1] = tvector3<T>(T(0), T(s), T(0));
			this->value[2] = tvector3<T>(T(0), T(0), T(s));
		}

		template <
			typename X1, typename Y1, typename Z1,
			typename X2, typename Y2, typename Z2,
			typename X3, typename Y3, typename Z3 >
			tmatrix3x3
			(
				X1 const & x1, Y1 const & y1, Z1 const & z1,
				X2 const & x2, Y2 const & y2, Z2 const & z2,
				X3 const & x3, Y3 const & y3, Z3 const & z3
			)
		{
			this->value[0] = col_type(T(x1), T(y1), T(z1));
			this->value[1] = col_type(T(x2), T(y2), T(z2));
			this->value[2] = col_type(T(x3), T(y3), T(z3));
		}

		template <typename V1, typename V2, typename V3>
		tmatrix3x3
		(
			tvector3<V1> const & v1,
			tvector3<V2> const & v2,
			tvector3<V3> const & v3
		)
		{
			this->value[0] = col_type(v1);
			this->value[1] = col_type(v2);
			this->value[2] = col_type(v3);
		}

		template <typename U>
		tmatrix3x3(tmatrix3x3<U> const & m)
		{
			this->value[0] = col_type(m[0]);
			this->value[1] = col_type(m[1]);
			this->value[2] = col_type(m[2]);
		}

		//-------------------- 运算符重载
		//-------------------- 求逆
		type _inverse() const
		{
			T S00 = value[0][0]; T S01 = value[0][1]; T S02 = value[0][2];

			T S10 = value[1][0]; T S11 = value[1][1]; T S12 = value[1][2];

			T S20 = value[2][0]; T S21 = value[2][1]; T S22 = value[2][2];

			//伴随矩阵
			type Inverse(
				S11 * S22 - S21 * S12, S12 * S20 - S22 * S10, S10 * S21 - S20 * S11,
				S02 * S21 - S01 * S22, S00 * S22 - S02 * S20, S01 * S20 - S00 * S21,
				S12 * S01 - S11 * S02, S10 * S02 - S12 * S00, S11 * S00 - S10 * S01);

			T Determinant = S00 * (S11 * S22 - S21 * S12) - S10 * (S01 * S22 - S21 * S02) + S20 * (S01 * S12 - S11 * S02);

			Inverse /= Determinant;
			return Inverse;
		}

		col_type & operator[](size_type i)
		{
			assert(i < this->length());
			return this->value[i];
		}

		col_type const & operator[](size_type i) const
		{
			assert(i < this->length());
			return this->value[i];
		}

		type & operator = (type const & m)
		{
			this->value[0] = m[0];
			this->value[1] = m[1];
			this->value[2] = m[2];
			return *this;
		}

		template <typename U>
		type & operator = (tmatrix3x3<U> const & m)
		{
			this->value[0] = m[0];
			this->value[1] = m[1];
			this->value[2] = m[2];
			return *this;
		}

		template <typename U>
		type & operator += (U const & s)
		{
			this->value[0] += s;
			this->value[1] += s;
			this->value[2] += s;
			return *this;
		}

		template <typename U>
		type & operator += (tmatrix3x3<U> const & m)
		{
			this->value[0] += m[0];
			this->value[1] += m[1];
			this->value[2] += m[2];
			return *this;
		}

		template <typename U>
		type & operator -= (U const & s)
		{
			this->value[0] -= s;
			this->value[1] -= s;
			this->value[2] -= s;
			return *this;
		}

		template <typename U>
		type & operator -= (tmatrix3x3<U> const & m)
		{
			this->value[0] -= m[0];
			this->value[1] -= m[1];
			this->value[2] -= m[2];
			return *this;
		}

		template <typename U>
		type & operator *= (U const & s)
		{
			this->value[0] *= s;
			this->value[1] *= s;
			this->value[2] *= s;
			return *this;
		}
		template <typename U>
		type & operator *= (tmatrix3x3<U> const & m)
		{
			return (*this = *this * m);
		}

		template <typename U>
		type & operator /= (U const & s)
		{
			this->value[0] /= s;
			this->value[1] /= s;
			this->value[2] /= s;
			return *this;
		}

		template <typename U>
		type & operator /= (tmatrix3x3<U> const & m)
		{
			return (*this = *this / m);
		}

		type & operator ++ ()
		{
			++this->value[0];
			++this->value[1];
			++this->value[2];
			return *this;
		}

		type & operator -- ()
		{
			--this->value[0];
			--this->value[1];
			--this->value[2];
			return *this;
		}

		tvector3<T> operator * (const tvector3<T> & v) const
		{
			return  tvector3<T>(
				value[0][0] * v[0] + value[1][0] * v[1] + value[2][0] * v[2],
				value[0][1] * v[0] + value[1][1] * v[1] + value[2][1] * v[2],
				value[0][2] * v[0] + value[1][2] * v[1] + value[2][2] * v[2]
				);
		}

		tvector2<T> operator * (const tvector2<T> & v) const
		{
			return  tvector2<T>(
				value[0][0] * v[0] + value[1][0] * v[1] + value[2][0],
				value[0][1] * v[0] + value[1][1] * v[1] + value[2][1]
				);
		}

		friend type operator + (type const & m, T const & s)
		{
			return type(
				m[0] + s,
				m[1] + s,
				m[2] + s);
		}

		friend type operator + (T const & s, type const & m)
		{
			return type(
				m[0] + s,
				m[1] + s,
				m[2] + s);
		}

		friend type operator + (type const & m1, type const & m2)
		{
			return type(
				m1[0] + m2[0],
				m1[1] + m2[1],
				m1[2] + m2[2]);
		}

		friend type operator - (type const & m, T const & s)
		{
			return type(
				m[0] - s,
				m[1] - s,
				m[2] - s);
		}

		friend type operator - (T const & s, type const & m)
		{
			return type(
				s - m[0],
				s - m[1],
				s - m[2]);
		}

		friend type operator - (type const & m1, type const & m2)
		{
			return type(
				m1[0] - m2[0],
				m1[1] - m2[1],
				m1[2] - m2[2]);
		}

		friend type operator * (type const & m, T const & s)
		{
			return type(
				m[0] * s,
				m[1] * s,
				m[2] * s);
		}

		friend type operator * (T const & s, type const & m)
		{
			return m * s;
		}

		friend tvector3<T> operator * (tvector3<T> const & v, type const & m)
		{
			return tvector3<T>(
				m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
				m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
				m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z);
		}

		friend type operator * (type const & m1, type const & m2)
		{
			T const srcA00 = m1[0][0]; T const srcA01 = m1[0][1]; T const srcA02 = m1[0][2];
			T const srcA10 = m1[1][0]; T const srcA11 = m1[1][1]; T const srcA12 = m1[1][2];
			T const srcA20 = m1[2][0]; T const srcA21 = m1[2][1]; T const srcA22 = m1[2][2];

			T const srcB00 = m2[0][0]; T const srcB01 = m2[0][1]; T const srcB02 = m2[0][2];
			T const srcB10 = m2[1][0]; T const srcB11 = m2[1][1]; T const srcB12 = m2[1][2];
			T const srcB20 = m2[2][0]; T const srcB21 = m2[2][1]; T const srcB22 = m2[2][2];

			type res;
			res[0][0] = srcA00 * srcB00 + srcA10 * srcB01 + srcA20 * srcB02;
			res[0][1] = srcA01 * srcB00 + srcA11 * srcB01 + srcA21 * srcB02;
			res[0][2] = srcA02 * srcB00 + srcA12 * srcB01 + srcA22 * srcB02;
			res[1][0] = srcA00 * srcB10 + srcA10 * srcB11 + srcA20 * srcB12;
			res[1][1] = srcA01 * srcB10 + srcA11 * srcB11 + srcA21 * srcB12;
			res[1][2] = srcA02 * srcB10 + srcA12 * srcB11 + srcA22 * srcB12;
			res[2][0] = srcA00 * srcB20 + srcA10 * srcB21 + srcA20 * srcB22;
			res[2][1] = srcA01 * srcB20 + srcA11 * srcB21 + srcA21 * srcB22;
			res[2][2] = srcA02 * srcB20 + srcA12 * srcB21 + srcA22 * srcB22;
			return res;
		}

		friend type operator / (type const & m, T const & s)
		{
			return type(
				m[0] / s,
				m[1] / s,
				m[2] / s);
		}

		friend type operator / (T const & s, type const & m)
		{
			return type(
				s / m[0],
				s / m[1],
				s / m[2]
			);
		}

		friend tvector3<T> operator / (type const & m, tvector3<T> const & v)
		{
			return m._inverse() * v;
		}

		friend tvector3<T> operator / (tvector3<T> const & v, type const & m)
		{
			return v * m._inverse();
		}

		friend type operator / (type const & m1, type const & m2)
		{
			return m1 * m2._inverse();
		}

		friend type const operator - (type const & m)
		{
			return type(
				-m[0],
				-m[1],
				-m[2]);
		}

		friend type const operator ++ (type const & m, int)
		{
			return type(
				m[0] + T(1),
				m[1] + T(1),
				m[2] + T(1));
		}

		friend type const operator -- (type const & m, int)
		{
			return type(
				m[0] - T(1),
				m[1] - T(1),
				m[2] - T(1));
		}

		friend bool operator == (type const & m1, type const & m2)
		{
			return (m1[0] == m2[0]) && (m1[1] == m2[1]) && (m1[2] == m2[2]);
		}

		friend bool operator != (type const & m1, type const & m2)
		{
			return (m1[0] != m2[0]) || (m1[1] != m2[1]) || (m1[2] != m2[2]);
		}

		//-------------------- 三种二维纹理基本操作
		//-------------------- 缩放
		void scale(T x, T y)
		{
			this->value[0] = col_type(T(x), T(0), T(0));
			this->value[1] = col_type(T(0), T(y), T(0));
			this->value[2] = col_type(T(0), T(0), T(1));
		}
		//-------------------- 绕坐标0,0旋转
		void rotate(T angle)
		{
			T rad = DEG2RAD(angle);
			T c = cos(rad);
			T s = sin(rad);
			this->value[0] = col_type(T(c), T(-s), T(0));
			this->value[1] = col_type(T(s), T(c), T(0));
			this->value[2] = col_type(T(0), T(0), T(1));
		}
		//-------------------- 平移
		void translate(T x, T y)
		{
			this->value[0] = col_type(T(1), T(0), T(0));
			this->value[1] = col_type(T(0), T(1), T(0));
			this->value[2] = col_type(T(x), T(y), T(1));
		}
	};

	//-------------------- 4x4 Matrix类（列优先）
	template <typename T>
	class tmatrix4x4
	{
	public:
		typedef T               value_type;
		typedef std::size_t     size_type;
		typedef tvector4<T>	col_type;
		typedef tvector4<T>	row_type;
		typedef tmatrix4x4<T>	type;
		typedef tmatrix4x4<T>	transpose_type;

	private:
		col_type value[4];

	public:
		size_type length() const
		{
			return 4;
		}
		size_type col_size()
		{
			return 4;
		}
		size_type row_size()
		{
			return 4;
		}

		void identify()
		{
			this->value[0] = col_type(1, 0, 0, 0);
			this->value[1] = col_type(0, 1, 0, 0);
			this->value[2] = col_type(0, 0, 1, 0);
			this->value[3] = col_type(0, 0, 0, 1);
		}

		col_type & operator[](size_type i)
		{
			assert(i < this->length());
			return this->value[i];
		}

		col_type const & operator[](size_type i) const
		{
			assert(i < this->length());
			return this->value[i];
		}

		T const * data() const
		{
			return &this->value[0][0];
		}

		//-------------------- 构造函数
		tmatrix4x4(type const & m)
		{
			this->value[0] = m.value[0];
			this->value[1] = m.value[1];
			this->value[2] = m.value[2];
			this->value[3] = m.value[3];
		}

		tmatrix4x4(tmatrix3x3<T> const & m)
		{
			this->value[0] = col_type(m[0], T(0));
			this->value[1] = col_type(m[1], T(0));
			this->value[2] = col_type(m[2], T(0));
			this->value[3] = col_type(T(0), T(0), T(0), T(1));
		}

		tmatrix4x4()
		{
			this->value[0] = col_type(T(1), T(0), T(0), T(0));
			this->value[1] = col_type(T(0), T(1), T(0), T(0));
			this->value[2] = col_type(T(0), T(0), T(1), T(0));
			this->value[3] = col_type(T(0), T(0), T(0), T(1));
		}

		tmatrix4x4(value_type s)
		{
			this->value[0] = col_type(s, T(0), T(0), T(0));
			this->value[1] = col_type(T(0), s, T(0), T(0));
			this->value[2] = col_type(T(0), T(0), s, T(0));
			this->value[3] = col_type(T(0), T(0), T(0), s);
		}

		tmatrix4x4
		(
			value_type const & x0, value_type const & y0, value_type const & z0, value_type const & w0,
			value_type const & x1, value_type const & y1, value_type const & z1, value_type const & w1,
			value_type const & x2, value_type const & y2, value_type const & z2, value_type const & w2,
			value_type const & x3, value_type const & y3, value_type const & z3, value_type const & w3
		)
		{
			this->value[0] = col_type(x0, y0, z0, w0);
			this->value[1] = col_type(x1, y1, z1, w1);
			this->value[2] = col_type(x2, y2, z2, w2);
			this->value[3] = col_type(x3, y3, z3, w3);
		}

		tmatrix4x4
		(
			col_type const & v0,
			col_type const & v1,
			col_type const & v2,
			col_type const & v3
		)
		{
			this->value[0] = v0;
			this->value[1] = v1;
			this->value[2] = v2;
			this->value[3] = v3;
		}

		template <typename U>
		tmatrix4x4(tmatrix4x4<U> const & m)
		{
			this->value[0] = col_type(m[0]);
			this->value[1] = col_type(m[1]);
			this->value[2] = col_type(m[2]);
			this->value[3] = col_type(m[3]);
		}

		template <typename U>
		tmatrix4x4(U const & s)
		{
			this->value[0] = tvector4<T>(T(s), T(0), T(0), T(0));
			this->value[1] = tvector4<T>(T(0), T(s), T(0), T(0));
			this->value[2] = tvector4<T>(T(0), T(0), T(s), T(0));
			this->value[3] = tvector4<T>(T(0), T(0), T(0), T(s));
		}
		template <
			typename X1, typename Y1, typename Z1, typename W1,
			typename X2, typename Y2, typename Z2, typename W2,
			typename X3, typename Y3, typename Z3, typename W3,
			typename X4, typename Y4, typename Z4, typename W4 >
			tmatrix4x4
			(
				X1 const & x1, Y1 const & y1, Z1 const & z1, W1 const & w1,
				X2 const & x2, Y2 const & y2, Z2 const & z2, W2 const & w2,
				X3 const & x3, Y3 const & y3, Z3 const & z3, W3 const & w3,
				X4 const & x4, Y4 const & y4, Z4 const & z4, W4 const & w4
			)
		{
			this->value[0] = col_type(T(x1), T(y1), T(z1), T(w1));
			this->value[1] = col_type(T(x2), T(y2), T(z2), T(w2));
			this->value[2] = col_type(T(x3), T(y3), T(z3), T(w3));
			this->value[3] = col_type(T(x4), T(y4), T(z4), T(w4));
		}

		template <typename V1, typename V2, typename V3, typename V4>
		tmatrix4x4
		(
			tvector4<V1> const & v1,
			tvector4<V2> const & v2,
			tvector4<V3> const & v3,
			tvector4<V4> const & v4
		)
		{
			this->value[0] = col_type(v1);
			this->value[1] = col_type(v2);
			this->value[2] = col_type(v3);
			this->value[3] = col_type(v4);
		}

		//-------------------- 运算符重载
		type & operator = (type const & m)
		{
			this->value[0] = m[0];
			this->value[1] = m[1];
			this->value[2] = m[2];
			this->value[3] = m[3];
			return *this;
		}

		template <typename U>
		type & operator = (tmatrix4x4<U> const & m)
		{
			this->value[0] = m[0];
			this->value[1] = m[1];
			this->value[2] = m[2];
			this->value[3] = m[3];
			return *this;
		}

		template <typename U>
		type & operator += (U const & s)
		{
			this->value[0] += s;
			this->value[1] += s;
			this->value[2] += s;
			this->value[3] += s;
			return *this;
		}

		template <typename U>
		type & operator += (tmatrix4x4<U> const & m)
		{
			this->value[0] += m[0];
			this->value[1] += m[1];
			this->value[2] += m[2];
			this->value[3] += m[3];
			return *this;
		}

		template <typename U>
		type & operator -= (U const & s)
		{
			this->value[0] -= s;
			this->value[1] -= s;
			this->value[2] -= s;
			this->value[3] -= s;
			return *this;
		}

		template <typename U>
		type & operator -= (tmatrix4x4<U> const & m)
		{
			this->value[0] -= m[0];
			this->value[1] -= m[1];
			this->value[2] -= m[2];
			this->value[3] -= m[3];
			return *this;
		}

		template <typename U>
		type & operator *= (U const & s)
		{
			this->value[0] *= s;
			this->value[1] *= s;
			this->value[2] *= s;
			this->value[3] *= s;
			return *this;
		}

		template <typename U>
		type & operator *= (tmatrix4x4<U> const & m)
		{
			return (*this = *this * m);
		}

		template <typename U>
		type & operator /= (U const & s)
		{
			this->value[0] /= s;
			this->value[1] /= s;
			this->value[2] /= s;
			this->value[3] /= s;
			return *this;
		}

		template <typename U>
		type & operator /= (tmatrix4x4<U> const & m)
		{
			return (*this = *this / m);
		}

		type & operator ++ ()
		{
			++this->value[0];
			++this->value[1];
			++this->value[2];
			++this->value[3];
			return *this;
		}

		type & operator -- ()
		{
			--this->value[0];
			--this->value[1];
			--this->value[2];
			--this->value[3];
			return *this;
		}

		friend tvector3<T> operator * (tvector3<T> const& v, tmatrix4x4<T> const& mat)
		{
			return tvector3<T>
				(
					v.x*mat[0][0] + v.y*mat[1][0] + v.z*mat[2][0] + 0 * mat[3][0],
					v.x*mat[0][1] + v.y*mat[1][1] + v.z*mat[2][1] + 0 * mat[3][1],
					v.x*mat[0][2] + v.y*mat[1][2] + v.z*mat[2][2] + 0 * mat[3][2]
					);
		}

		friend type operator + (type const & m, value_type s)
		{
			return type(
				m[0] + s,
				m[1] + s,
				m[2] + s,
				m[3] + s);
		}

		friend type operator + (value_type s, type const & m)
		{
			return type(
				m[0] + s,
				m[1] + s,
				m[2] + s,
				m[3] + s);
		}

		friend type operator + (type const & m1, type const & m2)
		{
			return type(
				m1[0] + m2[0],
				m1[1] + m2[1],
				m1[2] + m2[2],
				m1[3] + m2[3]);
		}

		friend type operator - (type const & m, value_type s)
		{
			return type(
				m[0] - s,
				m[1] - s,
				m[2] - s,
				m[3] - s);
		}

		friend type operator - (value_type s, type const & m)
		{
			return type(
				s - m[0],
				s - m[1],
				s - m[2],
				s - m[3]);
		}

		friend type operator - (type const & m1, type const & m2)
		{
			return type(
				m1[0] - m2[0],
				m1[1] - m2[1],
				m1[2] - m2[2],
				m1[3] - m2[3]);
		}

		friend type operator * (type const & m, value_type s)
		{
			return type(
				m[0] * s,
				m[1] * s,
				m[2] * s,
				m[3] * s);
		}

		friend type operator * (value_type s, type const & m)
		{
			return type(
				m[0] * s,
				m[1] * s,
				m[2] * s,
				m[3] * s);
		}

		friend col_type operator * (tmatrix4x4<T> const & m, row_type const & v)
		{
			return col_type(
				m[0][0] * v.x + m[1][0] * v.y + m[2][0] * v.z + m[3][0] * v.w,
				m[0][1] * v.x + m[1][1] * v.y + m[2][1] * v.z + m[3][1] * v.w,
				m[0][2] * v.x + m[1][2] * v.y + m[2][2] * v.z + m[3][2] * v.w,
				m[0][3] * v.x + m[1][3] * v.y + m[2][3] * v.z + m[3][3] * v.w);
		}

		friend row_type operator * (col_type const & v, tmatrix4x4<T> const & m)
		{
			return row_type(
				m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3] * v.w,
				m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3] * v.w,
				m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3] * v.w,
				m[3][0] * v.x + m[3][1] * v.y + m[3][2] * v.z + m[3][3] * v.w);
		}

		friend type operator * (type const & m1, type const & m2)
		{
			col_type const srcA0 = m1[0];
			col_type const srcA1 = m1[1];
			col_type const srcA2 = m1[2];
			col_type const srcA3 = m1[3];

			col_type const srcB0 = m2[0];
			col_type const srcB1 = m2[1];
			col_type const srcB2 = m2[2];
			col_type const srcB3 = m2[3];

			type res;
			res[0] = srcA0 * srcB0[0] + srcA1 * srcB0[1] + srcA2 * srcB0[2] + srcA3 * srcB0[3];
			res[1] = srcA0 * srcB1[0] + srcA1 * srcB1[1] + srcA2 * srcB1[2] + srcA3 * srcB1[3];
			res[2] = srcA0 * srcB2[0] + srcA1 * srcB2[1] + srcA2 * srcB2[2] + srcA3 * srcB2[3];
			res[3] = srcA0 * srcB3[0] + srcA1 * srcB3[1] + srcA2 * srcB3[2] + srcA3 * srcB3[3];
			return res;
		}

		friend type operator / (tmatrix4x4<T> const & m, value_type s)
		{
			return type(
				m[0] / s,
				m[1] / s,
				m[2] / s,
				m[3] / s);
		}

		friend type operator / (value_type s, tmatrix4x4<T> const & m)
		{
			return type(
				s / m[0],
				s / m[1],
				s / m[2],
				s / m[3]);
		}

		friend col_type operator / (tmatrix4x4<T> const & m, row_type const & v)
		{
			return m.inverse() * v;
		}

		friend row_type operator / (col_type const & v, tmatrix4x4<T> const & m)
		{
			return v * m.inverse();
		}

		friend tmatrix4x4<T> operator / (tmatrix4x4<T> const & m1, tmatrix4x4<T> const & m2)
		{
			return m1 * m2.inverse();
		}

		friend type const operator - (type const & m)
		{
			return type(
				-m[0],
				-m[1],
				-m[2],
				-m[3]);
		}

		friend type const operator ++ (type const & m, int)
		{
			return type(
				m[0] + T(1),
				m[1] + T(1),
				m[2] + T(1),
				m[3] + T(1));
		}

		friend type const operator -- (type const & m, int)
		{
			return type(
				m[0] - T(1),
				m[1] - T(1),
				m[2] - T(1),
				m[3] - T(1));
		}

		friend bool operator == (type const & m1, type const & m2)
		{
			return (m1[0] == m2[0]) && (m1[1] == m2[1]) && (m1[2] == m2[2]) && (m1[3] == m2[3]);
		}

		friend bool operator != (type const & m1, type const & m2)
		{
			return (m1[0] != m2[0]) || (m1[1] != m2[1]) || (m1[2] != m2[2]) || (m1[3] != m2[3]);
		}

		//-------------------- 矩阵操作
		//-------------------- 求逆
		type inverse() const
		{
			T subFactor00 = this->value[2][2] * this->value[3][3] - this->value[3][2] * this->value[2][3];
			T subFactor01 = this->value[2][1] * this->value[3][3] - this->value[3][1] * this->value[2][3];
			T subFactor02 = this->value[2][1] * this->value[3][2] - this->value[3][1] * this->value[2][2];
			T subFactor03 = this->value[2][0] * this->value[3][3] - this->value[3][0] * this->value[2][3];
			T subFactor04 = this->value[2][0] * this->value[3][2] - this->value[3][0] * this->value[2][2];
			T subFactor05 = this->value[2][0] * this->value[3][1] - this->value[3][0] * this->value[2][1];
			T subFactor06 = this->value[1][2] * this->value[3][3] - this->value[3][2] * this->value[1][3];
			T subFactor07 = this->value[1][1] * this->value[3][3] - this->value[3][1] * this->value[1][3];
			T subFactor08 = this->value[1][1] * this->value[3][2] - this->value[3][1] * this->value[1][2];
			T subFactor09 = this->value[1][0] * this->value[3][3] - this->value[3][0] * this->value[1][3];
			T subFactor10 = this->value[1][0] * this->value[3][2] - this->value[3][0] * this->value[1][2];
			T subFactor11 = this->value[1][1] * this->value[3][3] - this->value[3][1] * this->value[1][3];
			T SubFactor12 = this->value[1][0] * this->value[3][1] - this->value[3][0] * this->value[1][1];
			T subFactor13 = this->value[1][2] * this->value[2][3] - this->value[2][2] * this->value[1][3];
			T subFactor14 = this->value[1][1] * this->value[2][3] - this->value[2][1] * this->value[1][3];
			T subFactor15 = this->value[1][1] * this->value[2][2] - this->value[2][1] * this->value[1][2];
			T subFactor16 = this->value[1][0] * this->value[2][3] - this->value[2][0] * this->value[1][3];
			T subFactor17 = this->value[1][0] * this->value[2][2] - this->value[2][0] * this->value[1][2];
			T subFactor18 = this->value[1][0] * this->value[2][1] - this->value[2][0] * this->value[1][1];

			type res(
				+this->value[1][1] * subFactor00 - this->value[1][2] * subFactor01 + this->value[1][3] * subFactor02,
				-this->value[1][0] * subFactor00 + this->value[1][2] * subFactor03 - this->value[1][3] * subFactor04,
				+this->value[1][0] * subFactor01 - this->value[1][1] * subFactor03 + this->value[1][3] * subFactor05,
				-this->value[1][0] * subFactor02 + this->value[1][1] * subFactor04 - this->value[1][2] * subFactor05,

				-this->value[0][1] * subFactor00 + this->value[0][2] * subFactor01 - this->value[0][3] * subFactor02,
				+this->value[0][0] * subFactor00 - this->value[0][2] * subFactor03 + this->value[0][3] * subFactor04,
				-this->value[0][0] * subFactor01 + this->value[0][1] * subFactor03 - this->value[0][3] * subFactor05,
				+this->value[0][0] * subFactor02 - this->value[0][1] * subFactor04 + this->value[0][2] * subFactor05,

				+this->value[0][1] * subFactor06 - this->value[0][2] * subFactor07 + this->value[0][3] * subFactor08,
				-this->value[0][0] * subFactor06 + this->value[0][2] * subFactor09 - this->value[0][3] * subFactor10,
				+this->value[0][0] * subFactor11 - this->value[0][1] * subFactor09 + this->value[0][3] * SubFactor12,
				-this->value[0][0] * subFactor08 + this->value[0][1] * subFactor10 - this->value[0][2] * SubFactor12,

				-this->value[0][1] * subFactor13 + this->value[0][2] * subFactor14 - this->value[0][3] * subFactor15,
				+this->value[0][0] * subFactor13 - this->value[0][2] * subFactor16 + this->value[0][3] * subFactor17,
				-this->value[0][0] * subFactor14 + this->value[0][1] * subFactor16 - this->value[0][3] * subFactor18,
				+this->value[0][0] * subFactor15 - this->value[0][1] * subFactor17 + this->value[0][2] * subFactor18);

			T determinant =
				+this->value[0][0] * res[0][0]
				+ this->value[0][1] * res[1][0]
				+ this->value[0][2] * res[2][0]
				+ this->value[0][3] * res[3][0];

			res /= determinant;
			return res;
		}

		//-------------------- 平移
		type & translate(value_type x, value_type y, value_type z)
		{
			this->value[0] = col_type(1, 0, 0, 0);
			this->value[1] = col_type(0, 1, 0, 0);
			this->value[2] = col_type(0, 0, 1, 0);
			this->value[3] = col_type(x, y, z, 1);
			return  *this;
		}

		template<typename U>
		type & translate(U x, U y, U z)
		{
			this->value[0] = col_type(1, 0, 0, 0);
			this->value[1] = col_type(0, 1, 0, 0);
			this->value[2] = col_type(0, 0, 1, 0);
			this->value[3] = col_type(T(x), T(y), T(z), 1);
			return  *this;
		}

		type & translate(tvector3<T> const& pos)
		{
			this->value[0] = col_type(1, 0, 0, 0);
			this->value[1] = col_type(0, 1, 0, 0);
			this->value[2] = col_type(0, 0, 1, 0);
			this->value[3] = col_type(pos.x, pos.y, pos.z, 1);
			return  *this;
		}

		template<typename U>
		type & translate(tvector3<U> const& pos)
		{
			this->value[0] = col_type(1, 0, 0, 0);
			this->value[1] = col_type(0, 1, 0, 0);
			this->value[2] = col_type(0, 0, 1, 0);
			this->value[3] = col_type(T(pos.x), T(pos.y), T(pos.z), 1);
			return  *this;
		}

		//-------------------- 旋转
		type & rotate(value_type angle, tvector3<T> const & v)
		{
			T a = DEG2RAD(angle);
			T c = cos(a);
			T s = sin(a);

			tvector3<T> axis = normalize(v);

			tvector3<T> temp = (T(1) - c) * axis;

			type res;
			this->value[0][0] = c + temp[0] * axis[0];
			this->value[0][1] = 0 + temp[0] * axis[1] + s * axis[2];
			this->value[0][2] = 0 + temp[0] * axis[2] - s * axis[1];
			this->value[0][3] = 0;

			this->value[1][0] = 0 + temp[1] * axis[0] - s * axis[2];
			this->value[1][1] = c + temp[1] * axis[1];
			this->value[1][2] = 0 + temp[1] * axis[2] + s * axis[0];
			this->value[1][3] = 0;

			this->value[2][0] = 0 + temp[2] * axis[0] + s * axis[1];
			this->value[2][1] = 0 + temp[2] * axis[1] - s * axis[0];
			this->value[2][2] = c + temp[2] * axis[2];
			this->value[2][3] = 0;

			this->value[3][0] = 0;
			this->value[3][1] = 0;
			this->value[3][2] = 0;
			this->value[3][3] = 1;
			return  *this;
		}

		type & rotateX(value_type angle)
		{
			T a = DEG2RAD(angle);
			T c = cos(a);
			T s = sin(a);

			this->value[0] = col_type(1, 0, 0, 0);
			this->value[1] = col_type(0, c, s, 0);
			this->value[2] = col_type(0, -s, c, 0);
			this->value[3] = col_type(0, 0, 0, 1);

			return  *this;
		}

		template<typename U>
		type & rotateX(U angle)
		{
			T a = DEG2RAD(angle);
			T c = cos(a);
			T s = sin(a);

			this->value[0] = col_type(1, 0, 0, 0);
			this->value[1] = col_type(0, c, s, 0);
			this->value[2] = col_type(0, -s, c, 0);
			this->value[3] = col_type(0, 0, 0, 1);

			return  *this;
		}

		type & rotateY(value_type angle)
		{
			T a = DEG2RAD(angle);
			T c = cos(a);
			T s = sin(a);

			this->value[0] = col_type(c, 0, -s, 0);
			this->value[1] = col_type(0, 1, 0, 0);
			this->value[2] = col_type(s, 0, c, 0);
			this->value[3] = col_type(0, 0, 0, 1);
			return  *this;

		}

		template<typename U>
		type & rotateY(U angle)
		{
			T a = DEG2RAD(angle);
			T c = cos(a);
			T s = sin(a);

			this->value[0] = col_type(c, 0, -s, 0);
			this->value[1] = col_type(0, 1, 0, 0);
			this->value[2] = col_type(s, 0, c, 0);
			this->value[3] = col_type(0, 0, 0, 1);
			return  *this;

		}

		type & rotateZ(value_type angle)
		{
			T a = T(DEG2RAD(angle));
			T c = cos(a);
			T s = sin(a);

			this->value[0] = col_type(c, s, 0, 0);
			this->value[1] = col_type(-s, c, 0, 0);
			this->value[2] = col_type(0, 0, 1, 0);
			this->value[3] = col_type(0, 0, 0, 1);
			return  *this;
		}

		template<typename U>
		type & rotateZ(U angle)
		{
			T a = DEG2RAD(angle);

			T c = cos(a);
			T s = sin(a);

			this->value[0] = col_type(c, s, 0, 0);
			this->value[1] = col_type(-s, c, 0, 0);
			this->value[2] = col_type(0, 0, 1, 0);
			this->value[3] = col_type(0, 0, 0, 1);
			return  *this;
		}

		type rotateXY(T angleX, T angleY)
		{
			T cosX = cos(DEG2RAD(angleX));
			T sinX = sin(DEG2RAD(angleX));
			T cosY = cos(DEG2RAD(angleY));
			T sinY = sin(DEG2RAD(angleY));

			this->value[0] = col_type(cosY, -sinX * sinY, cosX * sinY, 0);
			this->value[1] = col_type(0, cosX, sinX, 0);
			this->value[2] = col_type(-sinY, -sinX * cosY, cosX * cosY, 0);
			this->value[3] = col_type(0, 0, 0, 1);
			return  *this;
		}

		type rotateYX(T angleX, T angleY)
		{
			T cosX = cos(DEG2RAD(angleX));
			T sinX = sin(DEG2RAD(angleX));
			T cosY = cos(DEG2RAD(angleY));
			T sinY = sin(DEG2RAD(angleY));

			this->value[0] = col_type(cosY, 0, sinY, 0);
			this->value[1] = col_type(-sinX * sinY, cosX, sinX * cosY, 0);
			this->value[2] = col_type(-cosX * sinY, -sinX, cosX * cosY, 0);
			this->value[3] = col_type(0, 0, 0, 1);

			return  *this;
		}

		type rotateYXZ(T yaw, T pitch, T roll)
		{
			T tmp_ch = cos(DEG2RAD(yaw));
			T tmp_sh = sin(DEG2RAD(yaw));
			T tmp_cp = cos(DEG2RAD(pitch));
			T tmp_sp = sin(DEG2RAD(pitch));
			T tmp_cb = cos(DEG2RAD(roll));
			T tmp_sb = sin(DEG2RAD(roll));

			type Result;
			this->value[0][0] = tmp_ch * tmp_cb + tmp_sh * tmp_sp * tmp_sb;
			this->value[0][1] = tmp_sb * tmp_cp;
			this->value[0][2] = -tmp_sh * tmp_cb + tmp_ch * tmp_sp * tmp_sb;
			this->value[0][3] = T(0);
			this->value[1][0] = -tmp_ch * tmp_sb + tmp_sh * tmp_sp * tmp_cb;
			this->value[1][1] = tmp_cb * tmp_cp;
			this->value[1][2] = tmp_sb * tmp_sh + tmp_ch * tmp_sp * tmp_cb;
			this->value[1][3] = T(0);
			this->value[2][0] = tmp_sh * tmp_cp;
			this->value[2][1] = -tmp_sp;
			this->value[2][2] = tmp_ch * tmp_cp;
			this->value[2][3] = T(0);
			this->value[3][0] = T(0);
			this->value[3][1] = T(0);
			this->value[3][2] = T(0);
			this->value[3][3] = T(1);

			return  *this;
		}

		/*
		*	pitch是围绕X轴旋转，也叫做俯仰角

		*	yaw是围绕Y轴旋转，也叫偏航角

		*	roll是围绕Z轴旋转，也叫翻滚角
		*/
		type yawPitchRoll(T yaw, T pitch, T roll)
		{
			T tmp_ch = cos(DEG2RAD(yaw));
			T tmp_sh = sin(DEG2RAD(yaw));
			T tmp_cp = cos(DEG2RAD(pitch));
			T tmp_sp = sin(DEG2RAD(pitch));
			T tmp_cb = cos(DEG2RAD(roll));
			T tmp_sb = sin(DEG2RAD(roll));

			this->value[0][0] = tmp_ch * tmp_cb + tmp_sh * tmp_sp * tmp_sb;
			this->value[0][1] = tmp_sb * tmp_cp;
			this->value[0][2] = -tmp_sh * tmp_cb + tmp_ch * tmp_sp * tmp_sb;
			this->value[0][3] = T(0);
			this->value[1][0] = -tmp_ch * tmp_sb + tmp_sh * tmp_sp * tmp_cb;
			this->value[1][1] = tmp_cb * tmp_cp;
			this->value[1][2] = tmp_sb * tmp_sh + tmp_ch * tmp_sp * tmp_cb;
			this->value[1][3] = T(0);
			this->value[2][0] = tmp_sh * tmp_cp;
			this->value[2][1] = -tmp_sp;
			this->value[2][2] = tmp_ch * tmp_cp;
			this->value[2][3] = T(0);
			this->value[3][0] = T(0);
			this->value[3][1] = T(0);
			this->value[3][2] = T(0);
			this->value[3][3] = T(1);

			return  *this;
		}

		//-------------------- 缩放
		type & scale(tvector3<T> const& s)
		{
			this->value[0] = col_type(s[0], T(0), T(0), T(0));
			this->value[1] = col_type(T(0), s[1], T(0), T(0));
			this->value[2] = col_type(T(0), T(0), s[2], T(0));
			this->value[3] = col_type(T(0), T(0), T(0), T(1));

			return  *this;
		}

		type & scale(value_type x, value_type y, value_type z)
		{
			this->value[0] = col_type(T(x), T(0), T(0), T(0));
			this->value[1] = col_type(T(0), T(y), T(0), T(0));
			this->value[2] = col_type(T(0), T(0), T(z), T(0));
			this->value[3] = col_type(T(0), T(0), T(0), T(1));

			return  *this;
		}

		template<typename U>
		type & scale(U x, U y, U z)
		{
			this->value[0] = col_type(T(x), T(0), T(0), T(0));
			this->value[1] = col_type(T(0), T(y), T(0), T(0));
			this->value[2] = col_type(T(0), T(0), T(z), T(0));
			this->value[3] = col_type(T(0), T(0), T(0), T(1));

			return  *this;
		}

		template<typename U, typename V, typename W>
		type & scale(U x, V y, W z)
		{
			this->value[0] = col_type(T(x), T(0), T(0), T(0));
			this->value[1] = col_type(T(0), T(y), T(0), T(0));
			this->value[2] = col_type(T(0), T(0), T(z), T(0));
			this->value[3] = col_type(T(0), T(0), T(0), T(1));

			return  *this;
		}

		type transpose() const
		{
			return  type(
				this->value[0][0], this->value[1][0], this->value[2][0], this->value[3][0],
				this->value[0][1], this->value[1][1], this->value[2][1], this->value[3][1],
				this->value[0][2], this->value[1][2], this->value[2][2], this->value[3][2],
				this->value[0][3], this->value[1][3], this->value[2][3], this->value[3][3]);
		}

		type extractMatrixRotation() const
		{
			return  type(
				this->value[0][0], this->value[0][1], this->value[0][2], 0.0,
				this->value[1][0], this->value[1][1], this->value[1][2], 0.0,
				this->value[2][0], this->value[2][1], this->value[2][2], 0.0,
				0.0, 0.0, 0.0, 1.0);
		}
	};

	//生成4x4旋转矩阵
	template <typename T>
	tmatrix4x4<T> rotateX(T angleX)
	{
		T cosX = cos(DEG2RAD(angleX));
		T sinX = sin(DEG2RAD(angleX));

		return tmatrix4x4<T>(
			T(1), T(0), T(0), T(0),
			T(0), cosX, sinX, T(0),
			T(0), -sinX, cosX, T(0),
			T(0), T(0), T(0), T(1));
	}

	template <typename T>
	tmatrix4x4<T> rotateY(T angleY)
	{
		T cosY = cos(DEG2RAD(angleY));
		T sinY = sin(DEG2RAD(angleY));

		return tmatrix4x4<T>(
			cosY, T(0), sinY, T(0),
			T(0), T(1), T(0), T(0),
			-sinY, T(0), cosY, T(0),
			T(0), T(0), T(0), T(1));
	}

	template <typename T>
	tmatrix4x4<T> rotateZ(T angleZ)
	{
		T cosZ = (T)cos(DEG2RAD(angleZ));
		T sinZ = (T)sin(DEG2RAD(angleZ));

		return tmatrix4x4<T>(
			cosZ, sinZ, T(0), T(0),
			-sinZ, cosZ, T(0), T(0),
			T(0), T(0), T(1), T(0),
			T(0), T(0), T(0), T(1));
	}

	template <typename T>
	tmatrix4x4<T> rotateXY(T angleX, T angleY)
	{
		return rotateX(angleX) * rotateY(angleY);
	}

	template <typename T>
	tmatrix4x4<T> rotateYX(T angleY, T angleX)
	{
		return rotateY(angleY) * rotateX(angleX);
	}

	template <typename T>
	tmatrix4x4<T> rotateXZ(T angleX, T angleZ)
	{
		return rotateX(angleX) * rotateZ(angleZ);
	}

	template <typename T>
	tmatrix4x4<T> rotateZX(T angleX, T angleZ)
	{
		return rotateZ(angleZ) * rotateX(angleX);
	}

	template <typename T>
	tmatrix4x4<T> rotateYXZ(T yaw, T pitch, T roll)
	{
		T tmp_ch = cos(DEG2RAD(yaw));
		T tmp_sh = sin(DEG2RAD(yaw));
		T tmp_cp = cos(DEG2RAD(pitch));
		T tmp_sp = sin(DEG2RAD(pitch));
		T tmp_cb = cos(DEG2RAD(roll));
		T tmp_sb = sin(DEG2RAD(roll));

		tmatrix4x4<T> res;
		res[0][0] = tmp_ch * tmp_cb + tmp_sh * tmp_sp * tmp_sb;
		res[0][1] = tmp_sb * tmp_cp;
		res[0][2] = -tmp_sh * tmp_cb + tmp_ch * tmp_sp * tmp_sb;
		res[0][3] = T(0);
		res[1][0] = -tmp_ch * tmp_sb + tmp_sh * tmp_sp * tmp_cb;
		res[1][1] = tmp_cb * tmp_cp;
		res[1][2] = tmp_sb * tmp_sh + tmp_ch * tmp_sp * tmp_cb;
		res[1][3] = T(0);
		res[2][0] = tmp_sh * tmp_cp;
		res[2][1] = -tmp_sp;
		res[2][2] = tmp_ch * tmp_cp;
		res[2][3] = T(0);
		res[3][0] = T(0);
		res[3][1] = T(0);
		res[3][2] = T(0);
		res[3][3] = T(1);
		return res;
	}

	template <typename T>
	tmatrix4x4<T> yawPitchRoll(T yaw, T pitch, T roll)
	{
		T tmp_ch = cos(DEG2RAD(yaw));
		T tmp_sh = sin(DEG2RAD(yaw));
		T tmp_cp = cos(DEG2RAD(pitch));
		T tmp_sp = sin(DEG2RAD(pitch));
		T tmp_cb = cos(DEG2RAD(roll));
		T tmp_sb = sin(DEG2RAD(roll));

		tmatrix4x4<T> res;
		res[0][0] = tmp_ch * tmp_cb + tmp_sh * tmp_sp * tmp_sb;
		res[0][1] = tmp_sb * tmp_cp;
		res[0][2] = -tmp_sh * tmp_cb + tmp_ch * tmp_sp * tmp_sb;
		res[0][3] = T(0);
		res[1][0] = -tmp_ch * tmp_sb + tmp_sh * tmp_sp * tmp_cb;
		res[1][1] = tmp_cb * tmp_cp;
		res[1][2] = tmp_sb * tmp_sh + tmp_ch * tmp_sp * tmp_cb;
		res[1][3] = T(0);
		res[2][0] = tmp_sh * tmp_cp;
		res[2][1] = -tmp_sp;
		res[2][2] = tmp_ch * tmp_cp;
		res[2][3] = T(0);
		res[3][0] = T(0);
		res[3][1] = T(0);
		res[3][2] = T(0);
		res[3][3] = T(1);
		return res;
	}

	//-------------------- 射线与三角形相交
	template<typename T>
	bool intersectTriangle(const tvector3<T>& orig, const tvector3<T>& dir, tvector3<T>& v0, tvector3<T>& v1, tvector3<T>& v2, T* t, T* u, T* v)
	{
		// Find vectors for two edges sharing vert0
		tvector3<T>    edge1 = v1 - v0;
		tvector3<T>    edge2 = v2 - v0;

		// Begin calculating determinant - also used to calculate U parameter
		tvector3<T>    pvec;
		pvec = cross(dir, edge2);

		// If determinant is near zero, ray lies in plane of triangle
		T   det = dot(edge1, pvec);

		tvector3<T>    tvec;
		if (det > 0)
		{
			tvec = orig - v0;
		}
		else
		{
			tvec = v0 - orig;
			det = -det;
		}
		if (det < 0.0001f)
			return false;
		// Calculate U parameter and test bounds
		*u = dot(tvec, pvec);
		if (*u < 0.0f || *u > det)
			return false;

		// Prepare to test V parameter
		tvector3<T>    qvec;
		qvec = cross(tvec, edge1);

		// Calculate V parameter and test bounds
		*v = dot(dir, qvec);
		if (*v < T(0.0f) || *u + *v > det)
			return false;

		*t = dot(edge2, qvec);
		T   fInvDet = T(1.0) / det;
		*t *= fInvDet;
		*u *= fInvDet;
		*v *= fInvDet;

		return true;
	}
	
	//-------------------- 计算三角形面积
	template<typename T> T calcTriangleArea(const tvector3<T>& pt1, const tvector3<T>& pt2, const tvector3<T>& pt3)
	{
		tvector3<T> e1 = pt2 - pt1;
		tvector3<T> e2 = pt3 - pt1;
		tvector3<T> e3 = cross(e1, e2);
		return  length(e3) * T(0.5);
	}

	//-------------------- 轴线角计算
	template <typename T>
	void axisAngle(tmatrix4x4<T> const & mat, tvector3<T> & axis, T & angle)
	{
		T epsilon = (T)0.01;
		T epsilon2 = (T)0.1;

		if ((fabs(mat[1][0] - mat[0][1]) < epsilon) &&
			(fabs(mat[2][0] - mat[0][2]) < epsilon) &&
			(fabs(mat[2][1] - mat[1][2]) < epsilon))
		{
			if ((fabs(mat[1][0] + mat[0][1]) < epsilon2) &&
				(fabs(mat[2][0] + mat[0][2]) < epsilon2) &&
				(fabs(mat[2][1] + mat[1][2]) < epsilon2) &&
				(fabs(mat[0][0] + mat[1][1] + mat[2][2] - (T)3.0) < epsilon2))
			{
				angle = (T)0.0;
				axis.x = (T)1.0;
				axis.y = (T)0.0;
				axis.z = (T)0.0;
				return;
			}
			angle = T(3.1415926535897932384626433832795);
			T xx = (mat[0][0] + (T)1.0) / (T)2.0;
			T yy = (mat[1][1] + (T)1.0) / (T)2.0;
			T zz = (mat[2][2] + (T)1.0) / (T)2.0;
			T xy = (mat[1][0] + mat[0][1]) / (T)4.0;
			T xz = (mat[2][0] + mat[0][2]) / (T)4.0;
			T yz = (mat[2][1] + mat[1][2]) / (T)4.0;
			if ((xx > yy) && (xx > zz))
			{
				if (xx < epsilon)
				{
					axis.x = (T)0.0;
					axis.y = (T)0.7071;
					axis.z = (T)0.7071;
				}
				else
				{
					axis.x = sqrt(xx);
					axis.y = xy / axis.x;
					axis.z = xz / axis.x;
				}
			}
			else if (yy > zz)
			{
				if (yy < epsilon)
				{
					axis.x = (T)0.7071;
					axis.y = (T)0.0;
					axis.z = (T)0.7071;
				}
				else
				{
					axis.y = sqrt(yy);
					axis.x = xy / axis.y;
					axis.z = yz / axis.y;
				}
			}
			else
			{
				if (zz < epsilon)
				{
					axis.x = (T)0.7071;
					axis.y = (T)0.7071;
					axis.z = (T)0.0;
				}
				else
				{
					axis.z = sqrt(zz);
					axis.x = xz / axis.z;
					axis.y = yz / axis.z;
				}
			}
			return;
		}
		T s = sqrt((mat[2][1] - mat[1][2]) * (mat[2][1] - mat[1][2]) + (mat[2][0] - mat[0][2]) * (mat[2][0] - mat[0][2]) + (mat[1][0] - mat[0][1]) * (mat[1][0] - mat[0][1]));
		if (abs(s) < T(0.001))
			s = (T)1.0;
		angle = acos((mat[0][0] + mat[1][1] + mat[2][2] - (T)1.0) / (T)2.0);
		axis.x = (mat[1][2] - mat[2][1]) / s;
		axis.y = (mat[2][0] - mat[0][2]) / s;
		axis.z = (mat[0][1] - mat[1][0]) / s;
	}

	//轴线角矩阵
	template <typename T>
	tmatrix4x4<T> axisAngleMatrix(tvector3<T> const & axis, T const angle)
	{
		T c = cos(angle);
		T s = sin(angle);
		T t = T(1) - c;
		tvector3<T> n = normalize(axis);

		return tmatrix4x4<T>(
			t * n.x * n.x + c, t * n.x * n.y + n.z * s, t * n.x * n.z - n.y * s, T(0),
			t * n.x * n.y - n.z * s, t * n.y * n.y + c, t * n.y * n.z + n.x * s, T(0),
			t * n.x * n.z + n.y * s, t * n.y * n.z - n.x * s, t * n.z * n.z + c, T(0),
			T(0), T(0), T(0), T(1));
	}

	//插值
	template <typename T>
	tmatrix4x4<T> interpolate(tmatrix4x4<T> const & m1, tmatrix4x4<T> const & m2, T const delta)
	{
		tmatrix4x4<T> m1rot = m1.extractMatrixRotation();
		tmatrix4x4<T> dltRotation = m2 * m1rot.transpose();
		dimension3<T> dltAxis;
		T dltAngle;
		axisAngle(dltRotation, dltAxis, dltAngle);
		tmatrix4x4<T> out = axisAngleMatrix(dltAxis, dltAngle * delta) * m1rot;
		out[3][0] = m1[3][0] + delta * (m2[3][0] - m1[3][0]);
		out[3][1] = m1[3][1] + delta * (m2[3][1] - m1[3][1]);
		out[3][2] = m1[3][2] + delta * (m2[3][2] - m1[3][2]);
		return out;
	}

	//正交投影
	template <typename valType>
	tmatrix4x4<valType> ortho(valType left, valType right, valType bottom, valType top, valType zNear, valType zFar)
	{
		tmatrix4x4<valType> res(1);
		res[0][0] = valType(2) / (right - left);
		res[1][1] = valType(2) / (top - bottom);
		res[2][2] = -valType(2) / (zFar - zNear);
		res[3][0] = -(right + left) / (right - left);
		res[3][1] = -(top + bottom) / (top - bottom);
		res[3][2] = -(zFar + zNear) / (zFar - zNear);
		return res;
	}

	//向量长度，距离，点乘，叉乘
	template <typename T>
	typename tvector2<T>::value_type length(tvector2<T> const & v)
	{
		typename tvector2<T>::value_type sqr = v.x * v.x + v.y * v.y;
		return sqrt(sqr);
	}

	template <typename T>
	typename tvector3<T>::value_type length(tvector3<T> const & v)
	{
		typename tvector3<T>::value_type sqr = v.x * v.x + v.y * v.y + v.z * v.z;
		return sqrt(sqr);
	}

	template <typename T>
	typename tvector4<T>::value_type length(tvector4<T> const & v)
	{
		typename tvector4<T>::value_type sqr = v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w;
		return sqrt(sqr);
	}

	template <typename T>
	typename tvector2<T>::value_type distance(tvector2<T> const & p0, tvector2<T> const & p1)
	{
		return length(p1 - p0);
	}

	template <typename T>
	typename tvector3<T>::value_type distance(tvector3<T> const & p0, tvector3<T> const & p1)
	{
		return length(p1 - p0);
	}

	template <typename T>
	typename tvector4<T>::value_type distance(tvector4<T> const & p0, tvector4<T> const & p1)
	{
		return length(p1 - p0);
	}

	template <typename T>
	typename tvector2<T>::value_type dot(tvector2<T> const & x, tvector2<T> const & y)
	{
		return x.x * y.x + x.y * y.y;
	}

	template <typename T>
	typename tvector3<T>::value_type dot(tvector3<T> const & x, tvector3<T> const & y)
	{
		return x.x * y.x + x.y * y.y + x.z * y.z;
	}

	template <typename T>
	typename tvector4<T>::value_type dot(tvector4<T> const & x, tvector4<T> const & y)
	{
		return x.x * y.x + x.y * y.y + x.z * y.z + x.w * y.w;
	}

	template <typename T>
	tvector3<T> cross(tvector3<T> const & x, tvector3<T> const & y)
	{
		return  tvector3<T>
			(
				x.y * y.z - y.y * x.z,
				x.z * y.x - y.z * x.x,
				x.x * y.y - y.x * x.y
				);
	}

	//向量单位化
	template <typename T>
	T inversesqrt(T x)
	{
		return T(1) / sqrt(x);
	}

	template <typename T>
	tvector2<T> normalize(tvector2<T> const & x)
	{
		typename tvector2<T>::value_type sqr = x.x * x.x + x.y * x.y;
		return tvector2<T>(x.x / sqr, x.y / sqr);
	}

	template <typename T>
	tvector3<T> normalize(tvector3<T> const & x)
	{
		typename tvector3<T>::value_type sqr = x.x * x.x + x.y * x.y + x.z * x.z;
		return x * inversesqrt(sqr);
	}

	template <typename T>
	tvector4<T> normalize(tvector4<T> const & x)
	{
		typename tvector4<T>::value_type sqr = x.x * x.x + x.y * x.y + x.z * x.z + x.w * x.w;
		return x * inversesqrt(sqr);
	}

	//生成投影矩阵 （行优先）（透视投影）
	template <typename valType>
	tmatrix4x4<valType> perspective(valType fovy, valType aspect, valType zNear, valType zFar)
	{
		valType range = tan(fovy * valType(DEG2RAD(0.5))) * zNear;
		valType left = -range * aspect;
		valType right = range * aspect;
		valType bottom = -range;
		valType top = range;

		tmatrix4x4<valType> res(valType(0));
		res[0][0] = (valType(2) * zNear) / (right - left);
		res[1][1] = (valType(2) * zNear) / (top - bottom);
		res[2][2] = -(zFar + zNear) / (zFar - zNear);
		res[2][3] = -valType(1);
		res[3][2] = -(valType(2) * zFar * zNear) / (zFar - zNear);
		return res;
	}

	//生成观察矩阵 （行优先）
	template <typename T>
	tmatrix4x4<T> lookAt(tvector3<T> const & eye, tvector3<T> const & center, tvector3<T> const & up) //人眼位置，向哪里观察，
	{
		tvector3<T> f = normalize(center - eye);
		tvector3<T> u = normalize(up);
		tvector3<T> s = normalize(cross(f, u));
		u = cross(s, f);

		tmatrix4x4<T> res(1);
		res[0][0] = s.x;
		res[1][0] = s.y;
		res[2][0] = s.z;
		res[0][1] = u.x;
		res[1][1] = u.y;
		res[2][1] = u.z;
		res[0][2] = -f.x;
		res[1][2] = -f.y;
		res[2][2] = -f.z;
		res[3][0] = -dot(s, eye);
		res[3][1] = -dot(u, eye);
		res[3][2] = dot(f, eye);
		return res;
	}

	//-------------------- 四元数类
	template<typename T>
	class tquat
	{
	public:
		typedef T value_type;
		typedef std::size_t size_type;

		value_type x;
		value_type y;
		value_type z;
		value_type w;

		size_type length() const
		{
			return	4;
		}

		tquat() :
			x(0),
			y(0),
			z(0),
			w(1) { }

		explicit tquat(value_type s, tvector3<T> const & v) :
			x(v.x),
			y(v.y),
			z(v.z),
			w(s)
		{
		}
		explicit tquat(tvector3<T> const & v, value_type s) :
			x(v.x),
			y(v.y),
			z(v.z),
			w(s)
		{
		}
		explicit tquat(value_type w, value_type x, value_type y, value_type z) :
			x(x),
			y(y),
			z(z),
			w(w)
		{}

		explicit tquat(tvector3<T> const & eulerAngle)
		{
			tvector3<T> c = cos(eulerAngle * value_type(0.5));
			tvector3<T> s = sin(eulerAngle * value_type(0.5));

			this->w = c.x * c.y * c.z + s.x * s.y * s.z;
			this->x = s.x * c.y * c.z - c.x * s.y * s.z;
			this->y = c.x * s.y * c.z + s.x * c.y * s.z;
			this->z = c.x * c.y * s.z - s.x * s.y * c.z;
		}

		explicit tquat(tmatrix3x3<T> const & m)
		{
			*this = quat_cast(m);
		}

		explicit tquat(tmatrix4x4<T> const & m)
		{
			*this = quat_cast(m);
		}

		value_type & operator[](int i)
		{
			return (&x)[i];
		}

		value_type const & operator[](int i) const
		{
			return (&x)[i];
		}

		tquat<T> & operator *= (value_type s)
		{
			this->w *= s;
			this->x *= s;
			this->y *= s;
			this->z *= s;
			return *this;
		}

		tquat<T> & operator = (const tquat<T>& right)
		{
			this->w = right.w;
			this->x = right.x;
			this->y = right.y;
			this->z = right.z;
			return *this;
		}

		tquat<T> & operator /= (value_type s)
		{
			this->w /= s;
			this->x /= s;
			this->y /= s;
			this->z /= s;
			return *this;
		}

		friend bool operator==(tquat<T> const & q1, tquat<T> const & q2)
		{
			return (q1.x == q2.x) && (q1.y == q2.y) && (q1.z == q2.z) && (q1.w == q2.w);
		}

		friend bool operator != (tquat<T> const & q1, tquat<T> const & q2)
		{
			return !(q1 == q2);
		}

		tquat<T> operator - ()
		{
			return tquat<T>(-w, -x, -y, -z);
		}

		friend tquat<T> operator+ (tquat<T> const & q, tquat<T> const & p)
		{
			return tquat<T>(
				q.w + p.w,
				q.x + p.x,
				q.y + p.y,
				q.z + p.z
				);
		}

		friend tquat<T> operator * (tquat<T> const & q, tquat<T> const & p)
		{
			return tquat<T>(
				q.w * p.w - q.x * p.x - q.y * p.y - q.z * p.z,
				q.w * p.x + q.x * p.w + q.y * p.z - q.z * p.y,
				q.w * p.y + q.y * p.w + q.z * p.x - q.x * p.z,
				q.w * p.z + q.z * p.w + q.x * p.y - q.y * p.x
				);
		}

		friend tvector3<T> operator * (tquat<T> const & q, tvector3<T> const & v)
		{
			value_type two(2);

			tvector3<T> uv;
			tvector3<T> uuv;
			tvector3<T> quatVector(q.x, q.y, q.z);
			uv = cross(quatVector, v);
			uuv = cross(quatVector, uv);
			uv *= two * q.w;
			uuv *= two;
			return v + uv + uuv;
		}


		friend tvector3<T> operator * (tvector3<T> const & v, tquat<T> const & q)
		{
			return  inverse(q) * v;
		}

		tquat<T> operator * (value_type s)
		{
			return tquat<T>(w * s, x * s, y * s, z * s);
		}

		friend tquat<T> operator * (value_type s, tquat<T> const & q)
		{
			return q * s;
		}

		tquat<T> operator / (value_type s)
		{
			return tquat<T>(w / s, x / s, y / s, z / s);
		}
	};

	template <typename T>
	T dot(tquat<T> const & q1, tquat<T> const & q2)
	{
		return q1.x * q2.x + q1.y * q2.y + q1.z * q2.z + q1.w * q2.w;
	}

	template <typename T>
	tquat<T> cross(tquat<T> const & q1, tquat<T> const & q2)
	{
		return tquat<T>(
			q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z,
			q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y,
			q1.w * q2.y + q1.y * q2.w + q1.z * q2.x - q1.x * q2.z,
			q1.w * q2.z + q1.z * q2.w + q1.x * q2.y - q1.y * q2.x);
	}

	template <typename T>
	T length(tquat<T> const & q)
	{
		return sqrt(dot(q, q));
	}

	template <typename T>
	tmatrix4x4<T> makeTransform(tvector3<T> const & position, tvector3<T> const& scale, const tquat<T>& orientation)
	{
		tmatrix3x3<T> rot3x3 = mat3_cast(orientation);

		return  tmatrix4x4<T>
			(
				scale.x * rot3x3[0][0], scale.x * rot3x3[0][1], scale.x * rot3x3[0][2], 0,
				scale.y * rot3x3[1][0], scale.y * rot3x3[1][1], scale.y * rot3x3[1][2], 0,
				scale.z * rot3x3[2][0], scale.z * rot3x3[2][1], scale.z * rot3x3[2][2], 0,
				position.x, position.y, position.z, 1
				);
	}

	template <typename T>
	T epsilon()
	{
		return std::numeric_limits<T>::epsilon();
	}

	template <typename T>
	tquat<T> inverse(tquat<T> const & q)
	{
		return conjugate(q) / dot(q, q);
	}

	template <typename T>
	tquat<T> conjugate(tquat<T> const & q)
	{
		return tquat<T>(q.w, -q.x, -q.y, -q.z);
	}

	template <typename T>
	tquat<T> mix(tquat<T> const & x, tquat<T> const & y, T const & a)
	{
		T cosTheta = dot(x, y);
		if (cosTheta > T(1) - epsilon<T>())
		{
			return tquat<T>(
				mix(x.w, y.w, a),
				mix(x.x, y.x, a),
				mix(x.y, y.y, a),
				mix(x.z, y.z, a)
				);
		}
		else
		{
			// Essential Mathematics, page 467
			T   angle = acos(cosTheta);
			return  (sin((T(1) - a) * angle) * x + sin(a * angle) * y) / sin(angle);
		}
	}

	template <typename T>
	tquat<T> lerp(tquat<T> const & x, tquat<T> const & y, T a)
	{
		assert(a >= T(0));
		assert(a <= T(1));
		return x * (T(1) - a) + (y * a);
	}

	template <typename T>
	tquat<T> slerp(tquat<T> const & x, tquat<T> const & y, T a)
	{
		tquat<T> z = y;

		T cosTheta = dot(x, y);

		if (cosTheta < T(0))
		{
			z = -y;
			cosTheta = -cosTheta;
		}
		if (cosTheta > T(1) - epsilon<T>())
		{
			return  tquat<T>
				(
					mix(x.w, z.w, a),
					mix(x.x, z.x, a),
					mix(x.y, z.y, a),
					mix(x.z, z.z, a)
					);
		}
		else
		{
			// Essential Mathematics, page 467
			T angle = acos(cosTheta);
			return (sin((T(1) - a) * angle) * x + sin(a * angle) * z) / sin(angle);
		}
	}

	template <typename T>
	tquat<T> rotate(typename tquat<T>::value_type angle, tvector3<T> const & axis)
	{
		tvector3<T> Tmp = axis;

		typename tquat<T>::value_type len = length(Tmp);
		if (abs(len - T(1)) > T(0.001f))
		{
			T oneOverLen = T(1) / len;
			Tmp.x *= oneOverLen;
			Tmp.y *= oneOverLen;
			Tmp.z *= oneOverLen;
		}
		typename tquat<T>::value_type const AngleRad = (T)DEG2RAD(angle);
		typename tquat<T>::value_type const Sin = (T)sin(AngleRad * T(0.5));
		return tquat<T>((T)cos(AngleRad * T(0.5)), Tmp.x * Sin, Tmp.y * Sin, Tmp.z * Sin);
	}

	template <typename valType>
	valType roll(tquat<valType> const & q)
	{
		return atan2(valType(2) * (q.x * q.y + q.w * q.z), q.w * q.w + q.x * q.x - q.y * q.y - q.z * q.z) * valType(RAD2DEG);
	}

	template <typename valType>
	valType pitch(tquat<valType> const & q)
	{
		return atan2(valType(2) * (q.y * q.z + q.w * q.x), q.w * q.w - q.x * q.x - q.y * q.y + q.z * q.z)* valType(RAD2DEG);
	}

	template <typename valType>
	valType yaw(tquat<valType> const & q)
	{
		return asin(valType(-2) * (q.x * q.z - q.w * q.y)) * valType(RAD2DEG);
	}

	template <typename T>
	tvector3<T> eulerAngles(tquat<T> const & x)
	{
		return tvector3<T>(pitch(x), yaw(x), roll(x));
	}

	template <typename T>
	tmatrix3x3<T> mat3_cast(const tquat<T>& q)
	{
		return tmatrix3x3<T>
			(
				1 - 2 * q.y * q.y - 2 * q.z * q.z, 2 * q.x * q.y + 2 * q.w * q.z, 2 * q.x * q.z - 2 * q.w * q.y,
				2 * q.x * q.y - 2 * q.w * q.z, 1 - 2 * q.x * q.x - 2 * q.z * q.z, 2 * q.y * q.z + 2 * q.w * q.x,
				2 * q.x * q.z + 2 * q.w * q.y, 2 * q.y * q.z - 2 * q.w * q.x, 1 - 2 * q.x * q.x - 2 * q.y * q.y
				);
	}

	template <typename T>
	tmatrix4x4<T> mat4_cast(tquat<T> const & q)
	{
		return tmatrix4x4<T>(mat3_cast(q));
	}

	template <typename T>
	tquat<T> quat_cast(tmatrix3x3<T> const & m)
	{
		typename tquat<T>::value_type fourXSquaredMinus1 = m[0][0] - m[1][1] - m[2][2];
		typename tquat<T>::value_type fourYSquaredMinus1 = m[1][1] - m[0][0] - m[2][2];
		typename tquat<T>::value_type fourZSquaredMinus1 = m[2][2] - m[0][0] - m[1][1];
		typename tquat<T>::value_type fourWSquaredMinus1 = m[0][0] + m[1][1] + m[2][2];

		int biggestIndex = 0;
		typename tquat<T>::value_type fourBiggestSquaredMinus1 = fourWSquaredMinus1;
		if (fourXSquaredMinus1 > fourBiggestSquaredMinus1)
		{
			fourBiggestSquaredMinus1 = fourXSquaredMinus1;
			biggestIndex = 1;
		}
		if (fourYSquaredMinus1 > fourBiggestSquaredMinus1)
		{
			fourBiggestSquaredMinus1 = fourYSquaredMinus1;
			biggestIndex = 2;
		}
		if (fourZSquaredMinus1 > fourBiggestSquaredMinus1)
		{
			fourBiggestSquaredMinus1 = fourZSquaredMinus1;
			biggestIndex = 3;
		}

		typename tquat<T>::value_type biggestVal = sqrt(fourBiggestSquaredMinus1 + T(1)) * T(0.5);
		typename tquat<T>::value_type mult = T(0.25) / biggestVal;

		tquat<T> res;
		switch (biggestIndex)
		{
		case 0:
			res.w = biggestVal;
			res.x = (m[1][2] - m[2][1]) * mult;
			res.y = (m[2][0] - m[0][2]) * mult;
			res.z = (m[0][1] - m[1][0]) * mult;
			break;
		case 1:
			res.w = (m[1][2] - m[2][1]) * mult;
			res.x = biggestVal;
			res.y = (m[0][1] + m[1][0]) * mult;
			res.z = (m[2][0] + m[0][2]) * mult;
			break;
		case 2:
			res.w = (m[2][0] - m[0][2]) * mult;
			res.x = (m[0][1] + m[1][0]) * mult;
			res.y = biggestVal;
			res.z = (m[1][2] + m[2][1]) * mult;
			break;
		case 3:
			res.w = (m[0][1] - m[1][0]) * mult;
			res.x = (m[2][0] + m[0][2]) * mult;
			res.y = (m[1][2] + m[2][1]) * mult;
			res.z = biggestVal;
			break;

		default:
			assert(false);
			break;
		}
		return res;
	}

	template <typename T>
	tquat<T> quat_cast(tmatrix4x4<T> const & m4)
	{
		return quat_cast(tmatrix3x3<T>(m4[0][0], m4[0][1], m4[0][2],
			m4[1][0], m4[1][1], m4[1][2],
			m4[2][0], m4[2][1], m4[2][2]));
	}

	template <typename T>
	T angle(tquat<T> const & x)
	{
		return acos(x.w) * T(2) * T(RAD2DEG);
	}

	template <typename T>
	tvector3<T> axis(tquat<T> const & x)
	{
		T   tmp1 = T(1) - x.w * x.w;
		if (tmp1 <= T(0))
		{
			return tvector3<T>(0, 0, 1);
		}
		T   tmp2 = T(1) / sqrt(tmp1);

		return tvector3<T>(x.x * tmp2, x.y * tmp2, x.z * tmp2);
	}

	template <typename valType>
	tquat<valType> angleAxis(valType angle, tvector3<valType> const & axis)
	{
		tquat<valType> result;

		valType a = (valType)(valType(DEG2RAD(angle)));
		valType s = sin(a * valType(0.5));

		result.w = cos(a * valType(0.5));
		result.x = axis.x * s;
		result.y = axis.y * s;
		result.z = axis.z * s;
		return result;
	}

	//-------------------- 轴向平行包围盒类3D
	template<typename T>
	class AxisAlignedBox
	{
	public:
		enum Extent
		{
			EXTENT_NULL,
			EXTENT_FINITE,
			EXTENT_INFINITE
		};

		tvector3<T>	_minimum;
		tvector3<T>	_maximum;
		Extent		_extent;

		/*
		   1-----2
		  /|    /|
		 / |   / |
		5-----4  |
		|  0--|--3
		| /   | /
		|/    |/
		6-----7
		*/
		typedef enum
		{
			FAR_LEFT_BOTTOM = 0,
			FAR_LEFT_TOP = 1,
			FAR_RIGHT_TOP = 2,
			FAR_RIGHT_BOTTOM = 3,
			NEAR_RIGHT_BOTTOM = 7,
			NEAR_LEFT_BOTTOM = 6,
			NEAR_LEFT_TOP = 5,
			NEAR_RIGHT_TOP = 4
		} CornerEnum;

		//-------------------- 构造函数
		AxisAlignedBox()
		{
			_minimum = tvector3<T>(T(-0.5), T(-0.5), T(-0.5));
			_maximum = tvector3<T>(T(0.5), T(0.5), T(0.5));
			_extent = EXTENT_NULL;
		}

		AxisAlignedBox(const AxisAlignedBox & rkBox)
		{
			setExtents(rkBox._minimum, rkBox._maximum);
			_extent = rkBox._extent;
		}

		AxisAlignedBox(const tvector3<T>& min, const tvector3<T>& max)
		{
			setExtents(min, max);
		}

		AxisAlignedBox(
			T mx, T my, T mz,
			T Mx, T My, T Mz
		)
		{
			setExtents(mx, my, mz, Mx, My, Mz);
		}

		AxisAlignedBox<T>& operator=(const AxisAlignedBox<T>& right)
		{
			setExtents(right._minimum, right._maximum);
			return *this;
		}

		//获取最小点坐标
		const tvector3<T>& getMinimum(void) const
		{
			return _minimum;
		}

		//获取最小点的可修改版本
		tvector3<T>& getMinimum(void)
		{
			return _minimum;
		}

		void setMinimum(const tvector3<T>& mins)
		{
			_minimum = mins;
		}
		void setMinimum(T x, T y, T z)
		{
			_minimum = tvector3<T>(x, y, z);
		}
		
		//获取最大点
		const tvector3<T>& getMaximum(void) const
		{
			return _maximum;
		}
		
		//获取最大值的可修改版本
		tvector3<T>& getMaximum(void)
		{
			return _maximum;
		}
		
		void setMaximum(const tvector3<T>& vec)
		{
			_maximum = vec;
		}

		void setMaximum(T x, T y, T z)
		{
			_maximum.x = x;
			_maximum.y = y;
			_maximum.z = z;
		}

		//调整框的一个尺寸的框的最大点的一个坐标
		void setMaximumX(T x)
		{
			_maximum.x = x;
		}

		void setMaximumY(T y)
		{
			_maximum.y = y;
		}

		void setMaximumZ(T z)
		{
			_maximum.z = z;
		}

		//一次性设置最大最小点
		void setExtents(const tvector3<T>& min, const tvector3<T>& max)
		{
			_minimum = min;
			_maximum = max;
			_extent = EXTENT_FINITE;
		}

		void setExtents(
			T mx, T my, T mz,
			T Mx, T My, T Mz)
		{
			_minimum.x = mx;
			_minimum.y = my;
			_minimum.z = mz;

			_maximum.x = Mx;
			_maximum.y = My;
			_maximum.z = Mz;
			_extent = EXTENT_FINITE;

		}

		/** 返回指向包含8个角点的数组的指针，对于有用碰撞与不对齐的对象。
		8个点顺序如下

		  1------2
		 /|    / |
		/ |   /  |
		5-----4  |
		|  0--|--3
		| /   | /
		|/    |/
		6-----7

		因为这个实现使用静态成员，一定要确保使用副本
		*/
		void getAllCorners(tvector3<T> mpCorners[8]) const
		{
			mpCorners[0] = _minimum;
			mpCorners[1].x = _minimum.x; mpCorners[1].y = _maximum.y; mpCorners[1].z = _minimum.z;
			mpCorners[2].x = _maximum.x; mpCorners[2].y = _maximum.y; mpCorners[2].z = _minimum.z;
			mpCorners[3].x = _maximum.x; mpCorners[3].y = _minimum.y; mpCorners[3].z = _minimum.z;

			mpCorners[4] = _maximum;
			mpCorners[5].x = _minimum.x; mpCorners[5].y = _maximum.y; mpCorners[5].z = _maximum.z;
			mpCorners[6].x = _minimum.x; mpCorners[6].y = _minimum.y; mpCorners[6].z = _maximum.z;
			mpCorners[7].x = _maximum.x; mpCorners[7].y = _minimum.y; mpCorners[7].z = _maximum.z;
		}

		//获取一个点的坐标
		tvector3<T> getCorner(CornerEnum cornerToGet) const
		{
			switch (cornerToGet)
			{
			case FAR_LEFT_BOTTOM:
				return _minimum;
			case FAR_LEFT_TOP:
				return tvector3<T>(_minimum.x, _maximum.y, _minimum.z);
			case FAR_RIGHT_TOP:
				return tvector3<T>(_maximum.x, _maximum.y, _minimum.z);
			case FAR_RIGHT_BOTTOM:
				return tvector3<T>(_maximum.x, _minimum.y, _minimum.z);
			case NEAR_RIGHT_BOTTOM:
				return tvector3<T>(_maximum.x, _minimum.y, _maximum.z);
			case NEAR_LEFT_BOTTOM:
				return tvector3<T>(_minimum.x, _minimum.y, _maximum.z);
			case NEAR_LEFT_TOP:
				return tvector3<T>(_minimum.x, _maximum.y, _maximum.z);
			case NEAR_RIGHT_TOP:
				return _maximum;
			default:
				return tvector3<T>();
			}
		}

		//将传入的框合并到当前框中，结果是包含两者的盒子
		void merge(const AxisAlignedBox<T>& right)
		{

			if ((right._extent == EXTENT_NULL) || (_extent == EXTENT_INFINITE))
			{
				return;
			}
			else if (right._extent == EXTENT_INFINITE)
			{
				_extent = EXTENT_INFINITE;
			}
			else if (_extent == EXTENT_NULL)
			{
				setExtents(right._minimum, right._maximum);
			}
			else
			{
				//! merge
				tvector3<T> min = _minimum;
				tvector3<T> max = _maximum;
				max.makeCeil(right._maximum);
				min.makeFloor(right._minimum);
				setExtents(min, max);
			}
		}

		//扩展框包含指定的点（如果需要）。
		void merge(const tvector3<T>& point)
		{
			switch (_extent)
			{
			case EXTENT_NULL: // if null, use this point
				setExtents(point, point);
				return;

			case EXTENT_FINITE:
				_maximum.makeCeil(point);
				_minimum.makeFloor(point);
				return;

			case EXTENT_INFINITE:
				return;
			}
		}

		void transform(const tmatrix4x4<T>& matrix)
		{
			tvector3<T>	oldMin;
			tvector3<T>	oldMax;
			tvector3<T>	currentCorner;

			oldMin = _minimum;
			oldMax = _maximum;


			//按照以下顺序依次计算角点：
			// 0，6，5，1，2，4，7，3
			//这个序列允许我们一次只更换一个成员以获得所有角落。

			//对于每一个，我们都使用矩阵进行转换
			//它给出结果点并合并结果点。

			currentCorner = oldMin;
			tvector3<T> vVert = currentCorner * matrix;
			setExtents(vVert, vVert);

			// First corner
			// min min min
			currentCorner = oldMin;
			merge(currentCorner * matrix);

			// min,min,max
			currentCorner.z = oldMax.z;
			merge(currentCorner * matrix);

			// min max max
			currentCorner.y = oldMax.y;
			merge(currentCorner * matrix);

			// min max min
			currentCorner.z = oldMin.z;
			merge(currentCorner * matrix);

			// max max min
			currentCorner.x = oldMax.x;
			merge(currentCorner * matrix);

			// max max max
			currentCorner.z = oldMax.z;
			merge(currentCorner * matrix);

			// max min max
			currentCorner.y = oldMin.y;
			merge(currentCorner * matrix);

			// max min min
			currentCorner.z = oldMin.z;
			merge(currentCorner * matrix);
		}

		//返回这个盒子是否与另一个盒子相交。
		bool intersects(const AxisAlignedBox& b2) const
		{
			if (_maximum.x < b2._minimum.x)
				return false;
			if (_maximum.y < b2._minimum.y)
				return false;
			if (_maximum.z < b2._minimum.z)
				return false;

			if (_minimum.x > b2._maximum.x)
				return false;
			if (_minimum.y > b2._maximum.y)
				return false;
			if (_minimum.z > b2._maximum.z)
				return false;
			return true;

		}

		//返回这个盒子是否与另一个盒子相交。不考虑Z
		bool intersectsNoZ(const AxisAlignedBox& b2) const
		{
			if (_maximum.x < b2._minimum.x)
				return false;
			if (_maximum.y < b2._minimum.y)
				return false;

			if (_minimum.x > b2._maximum.x)
				return false;
			if (_minimum.y > b2._maximum.y)
				return false;
			return true;

		}


		AxisAlignedBox<T> intersection(const AxisAlignedBox<T>& b2) const
		{
			tvector3<T> intMin = _minimum;
			tvector3<T> intMax = _maximum;

			intMin.makeCeil(b2.getMinimum());
			intMax.makeFloor(b2.getMaximum());

			if (intMin.x < intMax.x &&
				intMin.y < intMax.y &&
				intMin.z < intMax.z)
			{
				return AxisAlignedBox<T>(intMin, intMax);
			}

			return AxisAlignedBox<T>();
		}

		void setNull()
		{
			_extent = EXTENT_NULL;
		}

		bool isNull(void) const
		{
			return (_extent == EXTENT_NULL);
		}

		bool isFinite(void) const
		{
			return (_extent == EXTENT_FINITE);
		}

		void setInfinite()
		{
			_extent = EXTENT_INFINITE;
		}
		bool isInfinite(void) const
		{
			return (_extent == EXTENT_INFINITE);
		}

		void scale(const tvector3<T>& s)
		{
			tvector3<T> min = _minimum * s;
			tvector3<T> max = _maximum * s;
			setExtents(min, max);
		}

		bool intersects(const tvector3<T>& v) const
		{
			return(v.x >= _minimum.x  &&  v.x <= _maximum.x  &&
				v.y >= _minimum.y  &&  v.y <= _maximum.y  &&
				v.z >= _minimum.z  &&  v.z <= _maximum.z);
		}

		bool intersects(const tvector2<T>& v) const
		{
			return(v.x >= _minimum.x  &&  v.x <= _maximum.x  &&
				v.y >= _minimum.y  &&  v.y <= _maximum.y);
		}

		tvector3<T> getCenter(void) const
		{
			return tvector3<T>(
				(_maximum.x + _minimum.x) * T(0.5f),
				(_maximum.y + _minimum.y) * T(0.5f),
				(_maximum.z + _minimum.z) * T(0.5f)
				);
		}

		//获取大小
		tvector3<T> getSize(void) const
		{
			return _maximum - _minimum;
		}

		tvector3<T> getHalfSize(void) const
		{
			return (_maximum - _minimum) * T(0.5);
		}

		bool contains(const tvector3<T>& v) const
		{
			return _minimum.x <= v.x && v.x <= _maximum.x &&
				_minimum.y <= v.y && v.y <= _maximum.y &&
				_minimum.z <= v.z && v.z <= _maximum.z;
		}

		bool contains(const AxisAlignedBox& other) const
		{
			return this->_minimum.x <= other._minimum.x &&
				this->_minimum.y <= other._minimum.y &&
				this->_minimum.z <= other._minimum.z &&
				other._maximum.x <= this->_maximum.x &&
				other._maximum.y <= this->_maximum.y &&
				other._maximum.z <= this->_maximum.z;
		}
		bool operator== (const AxisAlignedBox& right) const
		{
			return this->_minimum == right._minimum &&
				this->_maximum == right._maximum;
		}
		bool operator!= (const AxisAlignedBox& right) const
		{
			return !(*this == right);
		}

	};

	//-------------------- 射线类
	template<typename T>
	class tray
	{
	private:
		tvector3<T>		_origin;
		tvector3<T>		_direction;

	public:
		typedef	T			value_type;
		typedef tray<T>		type;

		//-------------------- 构造函数
		tray() :
			_origin(T(0), T(0), T(0)),
			_direction(T(0), T(0), T(1)) { }

		tray(const tvector3<T>& origin, const tvector3<T>& direction) :
			_origin(origin),
			_direction(direction) { }


		void setOrigin(const tvector3<T>& origin)
		{
			_origin = origin;
		}

		const tvector3<T> & getOrigin(void) const
		{
			return _origin;
		}

		void setDirection(const tvector3<T>& dir)
		{
			_direction = dir;
		}

		const tvector3<T> & getDirection(void) const
		{
			return _direction;
		}

		//获取沿射线的点t单位的位置
		tvector3<T> getPoint(T time) const
		{
			return tvector3<T>(_origin + (_direction * time));
		}

		/**
		*   测试射线box相交
		*   如果相交,返回值中的first == true.否则false
		*   second为射线到点的距离
		*   调用getPoint方法，则返回交点
		*/
		std::pair<bool, T> intersects(const AxisAlignedBox<T>& box) const
		{
			T           lowt = 0.0f;
			T           t;
			bool        hit = false;
			tvector3<T>    hitpoint;
			tvector3<T>    min = box.getMinimum();
			tvector3<T>    max = box.getMaximum();

			/**
			*   点在包围盒里面
			*/
			if (_origin > min && _origin < max)
			{
				return std::pair<bool, T>(true, 0.0f);
			}

			// Check each face in turn, only check closest 3
			// Min x
			if (_origin.x <= min.x && _direction.x > 0)
			{
				t = (min.x - _origin.x) / _direction.x;
				if (t >= 0)
				{
					// Substitute t back into ray and check bounds and dist
					hitpoint = _origin + _direction * t;
					if (hitpoint.y >= min.y &&
						hitpoint.y <= max.y &&
						hitpoint.z >= min.z &&
						hitpoint.z <= max.z &&
						(!hit || t < lowt))
					{
						hit = true;
						lowt = t;
					}
				}
			}
			// Max x
			if (_origin.x >= max.x && _direction.x < 0)
			{
				t = (max.x - _origin.x) / _direction.x;
				if (t >= 0)
				{
					// Substitute t back into ray and check bounds and dist
					hitpoint = _origin + _direction * t;
					if (hitpoint.y >= min.y &&
						hitpoint.y <= max.y &&
						hitpoint.z >= min.z &&
						hitpoint.z <= max.z &&
						(!hit || t < lowt))
					{
						hit = true;
						lowt = t;
					}
				}
			}
			// Min y
			if (_origin.y <= min.y && _direction.y > 0)
			{
				t = (min.y - _origin.y) / _direction.y;
				if (t >= 0)
				{
					// Substitute t back into ray and check bounds and dist
					hitpoint = _origin + _direction * t;
					if (hitpoint.x >= min.x &&
						hitpoint.x <= max.x &&
						hitpoint.z >= min.z &&
						hitpoint.z <= max.z &&
						(!hit || t < lowt))
					{
						hit = true;
						lowt = t;
					}
				}
			}
			// Max y
			if (_origin.y >= max.y && _direction.y < 0)
			{
				t = (max.y - _origin.y) / _direction.y;
				if (t >= 0)
				{
					// Substitute t back into ray and check bounds and dist
					hitpoint = _origin + _direction * t;
					if (hitpoint.x >= min.x &&
						hitpoint.x <= max.x &&
						hitpoint.z >= min.z &&
						hitpoint.z <= max.z &&
						(!hit || t < lowt))
					{
						hit = true;
						lowt = t;
					}
				}
			}
			// Min z
			if (_origin.z <= min.z && _direction.z > 0)
			{
				t = (min.z - _origin.z) / _direction.z;
				if (t >= 0)
				{
					// Substitute t back into ray and check bounds and dist
					hitpoint = _origin + _direction * t;
					if (hitpoint.x >= min.x &&
						hitpoint.x <= max.x &&
						hitpoint.y >= min.y &&
						hitpoint.y <= max.y &&
						(!hit || t < lowt))
					{
						hit = true;
						lowt = t;
					}
				}
			}
			// Max z
			if (_origin.z >= max.z && _direction.z < 0)
			{
				t = (max.z - _origin.z) / _direction.z;
				if (t >= 0)
				{
					// Substitute t back into ray and check bounds and dist
					hitpoint = _origin + _direction * t;
					if (hitpoint.x >= min.x &&
						hitpoint.x <= max.x &&
						hitpoint.y >= min.y &&
						hitpoint.y <= max.y &&
						(!hit || t < lowt))
					{
						hit = true;
						lowt = t;
					}
				}
			}
			return std::pair<bool, T>(hit, lowt);
		}

	};

	//-------------------- 平面类（列优先）
	template <typename T>
	class tplane
	{
	public:
		tvector3<T>	_normal; //法向量
		T				_distance;

		//-------------------- 构造函数
		tplane()
		{
			_normal = tvector3<T>(0, 0, 0);
			_distance = 0.0f;
		}

		tplane(const tplane& right)
		{
			_normal = right._normal;
			_distance = right._distance;
		}

		//使用法向量以及沿法向量移动的距离
		tplane(const tvector3<T>& rkNormal, T fConstant)
		{
			_normal = rkNormal;
			_distance = -fConstant;
		}

		//使用4个数
		tplane(T x, T y, T z, T o)
		{
			_normal = tvector3<T>(x, y, z);
			T invLen = 1.0f / (_normal).length();
			_normal *= invLen;
			_distance = o * invLen;
		}

		//法向量及面上一点
		tplane(const tvector3<T>& rkNormal, const tvector3<T>& rkPoint)
		{
			redefine(rkNormal, rkPoint);
		}

		//面上三点
		tplane(const tvector3<T>& rkPoint0, const tvector3<T>& rkPoint1, const tvector3<T>& rkPoint2)
		{
			redefine(rkPoint0, rkPoint1, rkPoint2);
		}

		//-------------------- 
		//到点的距离
		float distance(const tvector3<T> &pos) const
		{
			return  dot(_normal, pos) + _distance;
		}

		bool operator==(const tplane & right) const
		{
			return (right._distance == _distance && right._normal == _normal);
		}
		bool operator!=(const tplane& right) const
		{
			return !(*this == right);
		}
	};

	//-------------------- 六个裁剪平面（截锥体）类 (列优先)
	template <typename T>
	class tfrustum
	{
	private:
		tplane<T> _planes[6];

	public:
		enum plane
		{
			FRUSTUM_LEFT = 0,
			FRUSTUM_RIGHT = 1,
			FRUSTUM_TOP = 2,
			FRUSTUM_BOTTOM = 3,
			FRUSTUM_FAR = 4,
			FRUSTUM_NEAR = 5,
		};

		//project * modle * view
		void loadFrustum(const tmatrix4x4<T> &mvp)
		{
			const T * dataPtr = mvp.data();
			_planes[FRUSTUM_LEFT] = tplane<T>(dataPtr[12] - dataPtr[0], dataPtr[13] - dataPtr[1], dataPtr[14] - dataPtr[2], dataPtr[15] - dataPtr[3]);
			_planes[FRUSTUM_RIGHT] = tplane<T>(dataPtr[12] + dataPtr[0], dataPtr[13] + dataPtr[1], dataPtr[14] + dataPtr[2], dataPtr[15] + dataPtr[3]);

			_planes[FRUSTUM_TOP] = tplane<T>(dataPtr[12] - dataPtr[4], dataPtr[13] - dataPtr[5], dataPtr[14] - dataPtr[6], dataPtr[15] - dataPtr[7]);
			_planes[FRUSTUM_BOTTOM] = tplane<T>(dataPtr[12] + dataPtr[4], dataPtr[13] + dataPtr[5], dataPtr[14] + dataPtr[6], dataPtr[15] + dataPtr[7]);

			_planes[FRUSTUM_FAR] = tplane<T>(dataPtr[12] - dataPtr[8], dataPtr[13] - dataPtr[9], dataPtr[14] - dataPtr[10], dataPtr[15] - dataPtr[11]);
			_planes[FRUSTUM_NEAR] = tplane<T>(dataPtr[12] + dataPtr[8], dataPtr[13] + dataPtr[9], dataPtr[14] + dataPtr[10], dataPtr[15] + dataPtr[11]);
		}

		//点
		bool pointInFrustum(const tvector3<T> &pos) const
		{
			for (int i = 0; i < 6; i++)
			{
				if (_planes[i].distance(pos) <= 0)
					return false;
			}
			return true;
		}

		//球
		bool sphereInFrustum(const tvector3<T> &pos, const float radius) const
		{
			for (int i = 0; i < 6; i++)
			{
				if (_planes[i].distance(pos) <= -radius)
					return false;
			}
			return true;
		}

		//立方体
		bool cubeInFrustum(T minX, T maxX, T minY, T maxY, T minZ, T maxZ) const
		{
			for (int i = 0; i < 6; i++)
			{
				if (_planes[i].distance(tvector3<T>(minX, minY, minZ)) > 0) continue;
				if (_planes[i].distance(tvector3<T>(minX, minY, maxZ)) > 0) continue;
				if (_planes[i].distance(tvector3<T>(minX, maxY, minZ)) > 0) continue;
				if (_planes[i].distance(tvector3<T>(minX, maxY, maxZ)) > 0) continue;
				if (_planes[i].distance(tvector3<T>(maxX, minY, minZ)) > 0) continue;
				if (_planes[i].distance(tvector3<T>(maxX, minY, maxZ)) > 0) continue;
				if (_planes[i].distance(tvector3<T>(maxX, maxY, minZ)) > 0) continue;
				if (_planes[i].distance(tvector3<T>(maxX, maxY, maxZ)) > 0) continue;
				return false;
			}
			return true;
		}

		const tplane<T> & getPlane(const int plane) const
		{
			return _planes[plane];
		}
	};

	//-------------------- RGBA颜色类
	class Rgba4Byte
	{
	public:
		Rgba4Byte(
			unsigned char r = 255,
			unsigned char g = 255,
			unsigned char b = 255,
			unsigned char a = 255
		)
			:_r(r), _g(g), _b(b), _a(a) { }
		Rgba4Byte(uint rgba)
		{
			_color = rgba;
		}
		bool operator == (const Rgba4Byte & right) const
		{
			return  _r == right._r &&
				_g == right._g &&
				_b == right._b &&
				_a == right._a;
		}
		bool operator != (const Rgba4Byte & right) const
		{
			return !(*this == right);
		}
		Rgba4Byte operator | (Rgba4Byte & right)
		{
			return Rgba4Byte(
				(_r + right._r) > 255 ? 255 : _r + right._r,
				(_g + right._g) > 255 ? 255 : _g + right._g,
				(_b + right._b) > 255 ? 255 : _b + right._b,
				(_a + right._a) > 255 ? 255 : _a + right._a
			);
		}
		friend Rgba4Byte operator + (const Rgba4Byte& left, const Rgba4Byte& right)
		{
			return Rgba4Byte(left._r * right._r
				, left._g * right._g
				, left._b * right._b
				, left._a * right._a);
		}
		operator unsigned()
		{
			return _color;
		}
		uint toUint()
		{
			return _color;
		}
		float length()
		{
			return sqrt(_r * _r + _g * _g + _b * _b) / (3.0 * 255 * 255);
		}
		union
		{
			struct
			{
				unsigned char _b;
				unsigned char _g;
				unsigned char _r;
				unsigned char _a;
			};
			uint _color;
		};

	};

	//-------------------- 颜色插值
	inline Rgba4Byte colorLerp(const Rgba4Byte& c1, const Rgba4Byte& c2, float s)
	{
		Rgba4Byte color;
		color._r = (unsigned char)(c1._r + s * (c2._r - c1._r));
		color._g = (unsigned char)(c1._g + s * (c2._g - c1._g));
		color._b = (unsigned char)(c1._b + s * (c2._b - c1._b));
		color._a = (unsigned char)(c1._a + s * (c2._a - c1._a));
		return color;
	}

	//-------------------- 坐标插值
	inline tvector2<float> uvLerp(const tvector2<float>& c1, const tvector2<float>& c2, float s)
	{
		tvector2<float> color;
		color.x = (c1.x + s * (c2.x - c1.x));
		color.y = (c1.y + s * (c2.y - c1.y));
		return  color;
	}

	typedef Rgba4Byte				Rgba;

	typedef tvector2<int>			int2;
	typedef tvector2<float>			float2;
	typedef tvector2<double>		double2;

	typedef tvector3<int>			int3;
	typedef tvector3<float>			float3;
	typedef tvector3<double>		double3;

	typedef tvector4<int>			int4;
	typedef tvector4<float>			float4;
	typedef tvector4<double>		double4;

	typedef tmatrix2x2<float>		matrix2;
	typedef tmatrix3x3<float>		matrix3;
	typedef tmatrix4x4<float>		matrix4;

	typedef tfrustum<float>			Frustum;
	typedef tray<float>				Ray;
	typedef tquat<float>			quatf;

	typedef AxisAlignedBox<float>	AABB3D;
}

#endif