#pragma once
#include <vector>
#include <type_traits>
#include <iostream>
#include <exception>
#include <algorithm>
#include <cmath>


namespace mtrx
{

	template <typename T>
	struct Matrix_cell_info;
	template <typename T>
	class Matrix;
	template <typename T>
	class Square_Matrix;

	template <typename T>
	void transform_matrix_with_main_element(mtrx::Matrix<T>&, std::vector<std::pair<size_t, T>>&, const size_t);
	template <typename T>
	size_t transform_square_matrix_with_main_element(mtrx::Square_Matrix<T>&, const size_t);
	
	template <typename T = double>
	struct Matrix_cell_info
	{
		T value;
		size_t x_pos{ 0 };
		size_t y_pos{ 0 };
		Matrix_cell_info() = default;
		~Matrix_cell_info() = default;
	};

	template <typename T = double>
	class Matrix final
	{
	private:
		size_t _size_x{ 0 };
		size_t _size_y{ 0 };
		std::vector<std::vector<T>> _data;
	public:
		Matrix(size_t x = 0, size_t y = 0) noexcept
			: _size_x(x), _size_y(y)
		{
			_data.resize(_size_y);
			for (auto& row : _data)
			{
				row.resize(_size_x);
				if (std::is_same<T, int>::value == false && std::is_same<T, double>::value == false && std::is_same<T, float>::value == false)
				{
					continue;
				}
				else
				{
					std::fill(row.begin(), row.end(), 0);
				}
			}
		}

		Matrix(const std::initializer_list<std::vector<T>>& list)
			: Matrix(list.begin()->size(), list.size())
		{
			size_t index = 0;
			for (auto& elem : list)
			{
				if (elem.size() != _size_x)
					throw std::exception("Error: length of rows is different");
				_data[index++] = elem;
			}
		}

		Matrix(const Matrix<T>& mat)
			: _size_x(mat._size_x), _size_y(mat._size_y)
		{
			_data.resize(_size_y);
			for (int i = 0; i < _size_y; i++)
			{
				_data[i] = mat._data[i];
			}
		}

		std::vector<T>& operator[](const size_t index) throw (std::out_of_range&)
		{
			if (index >= _size_y)
			{
				throw std::out_of_range("Error: out of range in Matrix");
			}
			return _data[index];
		}

		Matrix<T>& operator=(const Matrix<T>& new_m) noexcept
		{
			_size_x = new_m._size_x;
			_size_y = new_m._size_y;
			_data = new_m._data;
			return *this;
		}

		operator Square_Matrix<T>() const throw (std::exception&)
		{
			if (_size_x != _size_y)
			{
				throw std::exception("Error: it is not possible to convert a matrix to a square matrix");
			}
			Square_Matrix<T> sqr_mat(_size_x);
			for (int i = 0; i < _size_y; i++)
			{
				for (int j = 0; j < _size_x; j++)
				{
					sqr_mat[i][j] = get_xy(j, i);
				}
			}
			return sqr_mat;
		}

		Matrix_cell_info<T>&& find_max_abs_in(const size_t first_x, const size_t first_y, const size_t last_x, const size_t last_y) const throw (std::out_of_range&)
		{
			if (first_x >= _size_x || last_x >= _size_x || first_y >= _size_y || last_y >= _size_y)
			{
				throw std::out_of_range("Error: values for the search interval is out of range");
			}
			T max_value = abs(_data[first_y][first_x]);
			size_t X_pos = first_x;
			size_t y_pos = first_y;
			Matrix_cell_info<T> max_info;
			for (int i = first_y; i <= last_y; i++)
			{
				for (int j = first_x; j <= last_x; j++)
				{
					if (abs(_data[i][j]) >= max_value)
					{
						max_info.value = abs(_data[i][j]);
						max_info.x_pos = j;
						max_info.y_pos = i;
						max_value = max_info.value;
					}
				}
			}
			return std::move(max_info);
		}

		Matrix_cell_info<T>&& find_max_abs() const noexcept
		{
			return std::forward(find_max_abs_in(0, 0, _size_x, _size_y));
		}

		Matrix<T>& swap_rows(const size_t ind_1, const size_t ind_2) throw (std::out_of_range&)
		{
			std::swap(_data[ind_1], _data[ind_2]);
			return *this;
		}

		Matrix<T>& swap_columns(const size_t ind_1, const size_t ind_2) throw (std::out_of_range&)
		{
			std::for_each(_data.begin(), _data.end(), [=](auto& row){ std::swap(row.at(ind_1), row.at(ind_2)); });
			return *this;
		}

		size_t row_size() const noexcept
		{
			return _size_x;
		}

		size_t column_size() const noexcept
		{
			return _size_y;
		}

		void resize_x(const size_t new_size_x) throw (std::exception&)
		{
			if (new_size_x < _size_x)
			{
				throw std::exception("Error: new size in resize_x is less than old size");
			}
			for (auto& row : _data)
			{
				row.resize(new_size_x);
				if (std::is_same<T, int>::value == false && std::is_same<T, double>::value == false && std::is_same<T, float>::value == false)
				{
					continue;
				}
				else
				{
					std::fill(row.begin() + new_size_x, row.end(), 0);
				}
			}
			_size_x = new_size_x;
		}

		void resize_y(const size_t new_size_y) throw (std::exception&)
		{
			if (new_size_y < _size_y)
			{
				throw std::exception("Error: new size in resize_y is less than old size");
			}
			_data.resize(new_size_y);
			for (int i = _size_y; i < new_size_y; i++)
			{
				_data[i].resize(_size_x);
				if (std::is_same<T, int>::value == false && std::is_same<T, double>::value == false && std::is_same<T, float>::value == false)
				{
					continue;
				}
				else
				{
					std::fill(_data[i].begin(), _data[i].end(), 0);
				}
			}
			_size_y = new_size_y;
		}

		T get_xy(const size_t x, const size_t y) const throw (std::out_of_range&)
		{
			if (x >= _size_x || y >= _size_y)
			{
				throw std::out_of_range("Error: out of range in get_xy");
			}
			return _data[y][x];
		}

		Matrix<T> operator-() const noexcept
		{
			Matrix<T> new_m = *this;
			for (auto& row : new_m._data)
			{
				for (auto& elem : row)
				{
					elem = -elem;
				}
			}
			return new_m;
		}

		~Matrix() = default;
		
		friend std::ostream& operator<<(std::ostream& out, const Matrix<T>& m)
		{
			for (const auto& row : m._data)
			{
				for (const auto& elem : row)
				{
					out.width(12);
					out << elem;
				}
				out << std::endl;
			}
			return out;
		}
	};

	template <typename T = double>
	class Square_Matrix final
	{
	private:
		Matrix<T> _matrix;

		void get_matr_without_ij(Square_Matrix& matr, int i, int j) const
		{
			for (int k = 0; k < size() - 1; k++)
			{
				for (int r = 0; r < size() - 1; r++)
				{
					if (r < j && k < i)
						matr._matrix[k][r] = _matrix.get_xy(r, k);
					if (r >= j && k < i)
						matr._matrix[k][r] = _matrix.get_xy(r + 1, k);
					if (r < j && k >= i)
						matr._matrix[k][r] = _matrix.get_xy(r, k + 1);
					if (r >= j && k >= i)
						matr._matrix[k][r] = _matrix.get_xy(r + 1, k + 1);
				}
			}
		}

	public:
		Square_Matrix(size_t size = 0) noexcept
			: _matrix(size, size)
		{}

		Square_Matrix(const std::initializer_list<std::vector<T>>& list) noexcept
			: _matrix(list)
		{}

		Square_Matrix(const Square_Matrix<T>& s_mat) noexcept
			: _matrix(s_mat._matrix)
		{}

		operator Matrix<T>() const
		{
			return _matrix;
		}

		Square_Matrix<T>& operator=(const Square_Matrix<T>& new_m) noexcept
		{
			_matrix = new_m._matrix;
			return *this;
		}

		std::vector<T>& operator[](const size_t index) throw (std::out_of_range&)
		{
			return _matrix[index];
		}

		Square_Matrix<T> operator-() const noexcept
		{
			auto new_m = *this;
			new_m._matrix = -new_m._matrix;
			return new_m;
		}

		Matrix_cell_info<T> find_max_abs_in(const size_t first_x, const size_t first_y, const size_t last_x, const size_t last_y) const throw (std::out_of_range&)
		{
			return _matrix.find_max_abs_in(first_x, first_y, last_x, last_y);
		}

		Matrix_cell_info<T> find_max_abs() const noexcept
		{
			return _matrix.find_max_abs_in(0, 0, _matrix.row_size(), _matrix.column_size());
		}

		Square_Matrix<T>& swap_rows(const size_t ind_1, const size_t ind_2) throw (std::out_of_range&)
		{
			_matrix.swap_rows(ind_1, ind_2);
			return *this;
		}

		Square_Matrix<T>& swap_columns(const size_t ind_1, const size_t ind_2) throw (std::out_of_range&)
		{
			_matrix.swap_columns(ind_1, ind_2);
			return *this;
		}

		size_t size() const noexcept
		{
			return _matrix.column_size();
		}

		Matrix<T> get_Matrix() const noexcept
		{
			return _matrix;
		}

		double det() const noexcept
		{
			Square_Matrix<T> copy = *this;
			//std::cout << copy << std::endl;
			size_t count_of_swap = 0;
			for (int i = 0; i < copy.size(); i++)
			{
				count_of_swap += transform_square_matrix_with_main_element(copy, i);
				//std::cout << copy << std::endl;
				for (int j = i + 1; j < copy.size(); j++)
				{
					if (std::abs(copy[i][i]) < 10e-10)
					{
						return 0;
					}
					double q = static_cast<double>(copy[j][i]) / copy[i][i];
					for (int k = i; k < copy.size(); k++)
					{
						copy[j][k] -= copy[i][k] * q;
					}
				}
				//std::cout << copy << std::endl;
			}
			double det = pow(-1, count_of_swap);
			
			for (int i = 0; i < copy.size(); i++)
			{
				det *= copy[i][i];
			}
			return det;
		}

		Square_Matrix<T> inverse() const throw (std::exception&)
		{
			//std::cout << "-----------------------------------------------" << std::endl;
			//std::cout << "Matrix:" << std::endl << *this << std::endl;
			Square_Matrix<T> inv(size());
			Square_Matrix<T> minor(size() - 1);
			double determ = det() + 1;      //  т.к. double нельзя сравнивать с 0
			if (determ == 1)
			{
				throw std::exception("Error: determinant is equal to 0");
			}
			determ--;
			for (int i = 0, sign_i = 1; i < size(); i++, sign_i = -sign_i)
			{
				for (int j = 0, sign_j = 1; j < size(); j++, sign_j = -sign_j)
				{
					get_matr_without_ij(minor, i, j);
					inv._matrix[j][i] = static_cast<double>(sign_i) * static_cast<double>(sign_j) * minor.det() / determ;
					//std::cout << "i = " << i << "\t" << "j = " << j << std::endl << "minor = " << std::endl << minor << std::endl << "det = " << minor.det() << std::endl;
					//std::cout << "elem[j][i] = " << inv._matrix[j][i] << std::endl;
				}
			}
			//std::cout << "------------------------------------------------" << std::endl;
			return inv;
		}

		T get_xy(const size_t x, const size_t y) const throw (std::exception&)
		{
			return _matrix.get_xy(x, y);
		}

		~Square_Matrix() = default;
		
		friend std::ostream& operator<<(std::ostream& out, const Square_Matrix<T>& m)
		{
			out << m._matrix;
			return out;
		}
	};

	template <typename T = double>
	Matrix<T> make_expanded_matrix(const Square_Matrix<T>& sq_mat, const std::vector<T>& f)
	{
		Matrix<T> exp_matrix = sq_mat.get_Matrix();
		for (int i = 0; i < exp_matrix.column_size(); i++)
		{
			exp_matrix.resize_x(sq_mat.size() + 1);
			exp_matrix[i][exp_matrix.row_size() - 1] = f[i];
		}
		return exp_matrix;
	}

	template<typename T>
	std::vector<T> operator*(const Matrix<T>& mat, const std::vector<T>& vec) throw (std::exception&)
	{
		if (mat.row_size() != vec.size())
		{
			throw std::exception("Error: incorrect dimensions of the elements being multiplied");
		}
		std::vector<T> result(mat.column_size());
		for (int i = 0; i < result.size(); i++)
		{
			result[i] = 0;
			for (int j = 0; j < vec.size(); j++)
			{
				result[i] += mat.get_xy(j, i) * vec[j];
			}
		}
		return result;
	}

	template<typename T>
	std::vector<T> operator*(const Square_Matrix<T>& mat, const std::vector<T>& vec) throw (std::exception&)
	{
		return mat.get_Matrix() * vec;
	}

	template<typename T>
	Matrix<T> operator*(const Matrix<T>& mat_1, const Matrix<T> mat_2) throw (std::exception&)
	{
		if (mat_1.row_size() != mat_2.column_size())
		{
			throw std::exception("Error: incorrect dimensions of the elements being multiplied");
		}
		Matrix<T> result(mat_1.column_size(), mat_2.row_size());
		for (int i = 0; i < mat_1.column_size(); i++)
		{
			for (int j = 0; j < mat_2.row_size(); j++)
			{
				result[i][j] = mat_1.get_xy(0, i) * mat_2.get_xy(j, 0);
				for (int k = 1; k < mat_1.row_size(); k++)
				{
					result[i][j] += mat_1.get_xy(k, i) * mat_2.get_xy(j, k);
				}
			}
		}
		return result;
	}

	template <typename T>
	Square_Matrix<T> operator*(const Square_Matrix<T>& m_1, const Square_Matrix<T>& m_2) throw (std::exception&)
	{
		return m_1.get_Matrix() * m_2.get_Matrix();
	}

	template <typename T>
	Matrix<T> operator+(const Matrix<T>& mat_1, const Matrix<T> mat_2) throw (std::exception&)
	{
		if ((mat_1.row_size() != mat_1.row_size()) || (mat_1.column_size() != mat_2.column_size()))
		{
			throw std::exception("Error: incorrect dimensions of the elements being summarized");
		}
		Matrix<T> result(mat_1.column_size(), mat_1.row_size());
		for (int i = 0; i < mat_1.column_size(); i++)
		{
			for (int j = 0; j < mat_1.row_size(); j++)
			{
				result[i][j] = mat_1.get_xy(j, i) + mat_2.get_xy(j, i);
			}
		}
		return result;
	}

	template <typename T>
	Square_Matrix<T> operator+(const Square_Matrix<T>& m_1, const Square_Matrix<T>& m_2) throw (std::exception&)
	{
		return m_1.get_Matrix() + m_2.get_Matrix();
	}

	template <typename T>
	double norm(std::vector<T>& v) noexcept
	{
		double sum = 0;
		for (auto& elem : v)
		{
			sum += elem * elem;
		}
		return sqrt(sum);
	}

	//	Преобразование матрицы (перемещение максимального по модулю элемента в нужную позицию)

	template <typename T = double>
	void transform_matrix_with_main_element(mtrx::Matrix<T>& matrix, std::vector<std::pair<size_t, T>>& var, const size_t num_first)
	{
		auto info_about_max = matrix.find_max_abs_in(num_first, num_first, matrix.row_size() - 2, matrix.column_size() - 1);
		if (info_about_max.y_pos != num_first)
			matrix.swap_rows(info_about_max.y_pos, num_first);
		if (info_about_max.x_pos != num_first)
		{
			matrix.swap_columns(info_about_max.x_pos, num_first);
			std::swap(var[info_about_max.x_pos].first, var[num_first].first);
		}
	}

	template <typename T>
	size_t transform_square_matrix_with_main_element(mtrx::Square_Matrix<T>& mat, const size_t num_first)
	{
		auto info_about_max = mat.find_max_abs_in(num_first, num_first, mat.size() - 1, mat.size() - 1);
		size_t count_of_swap = 0;
		if (info_about_max.y_pos != num_first)
		{
			mat.swap_rows(info_about_max.y_pos, num_first);
			count_of_swap++;
		}
		if (info_about_max.x_pos != num_first)
		{
			mat.swap_columns(info_about_max.x_pos, num_first);
			count_of_swap++;
		}
		return count_of_swap;
	}


}	// namespace mtrx

template <typename T>
std::vector<T> operator+(const std::vector<T>& v_1, const std::vector<T>& v_2) throw (std::exception&)
{
	if (v_1.size() != v_2.size())
	{
		throw std::exception("Error: incorrect dimensions of the elements being subtracted");
	}
	auto result = v_1;
	for (int i = 0; i < result.size(); i++)
	{
		result[i] += v_2[i];
	}
	return result;
}

template <typename T>
std::vector<T> operator-(const std::vector<T>& v_1, const std::vector<T>& v_2) throw (std::exception&)
{
	if (v_1.size() != v_2.size())
	{
		throw std::exception("Error: incorrect dimensions of the elements being subtracted");
	}
	auto result = v_1;
	for (int i = 0; i < result.size(); i++)
	{
		result[i] -= v_2[i];
	}
	return result;
}