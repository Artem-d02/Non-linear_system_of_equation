#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <exception>
#include <fstream>

#include "Matrix.h"

//-------------------------------- class for information about method ---------------

struct info_about_method
{
	std::vector<double> solution;
	size_t number_of_it{ 0 };
};

std::vector<double> system_operator(const std::vector<double>&) noexcept;

info_about_method DIM(std::vector<double>, double, std::function<std::vector<double>(const std::vector<double>&)>) noexcept;

info_about_method Newton_method(std::vector<double>, double, std::function<mtrx::Square_Matrix<double>(const std::vector<double>&)>, std::function<std::vector<double>(const std::vector<double>&)>) noexcept;

template <typename T>
std::vector<T> operator-(const std::vector<T>&, const std::vector<T>&);

template <typename T>
std::ostream& operator<<(std::ostream&, const std::vector<T>&);

std::vector<double> F_of_system(const std::vector<double>&) noexcept;

mtrx::Square_Matrix<double> W_of_system(const std::vector<double>&) noexcept;

//-------------------------------- main ---------------------------------------------
// 
//-------------------------------- System: ------------------------------------------
//------------------------------ { sin(x) - y = 1.32 --------------------------------
//------------------------------ { cos(y) - x = -0.85 -------------------------------

int main()
{
	std::fstream fout;
	fout.open("Output.txt", std::fstream::out);	

	std::vector<double> start = { 0.0, 0.0 };		// initial approximation
	double E = 1e-5;								// accuracy

	std::cout << "Initial approximation:" << std::endl << start;
	std::cout << "Accuracy of solution: " << E << std::endl << std::endl;

	fout << "Initial approximation:" << std::endl << start;
	fout << "Accuracy of solution: " << E << std::endl << std::endl;

	auto sol_info = DIM(start, E, system_operator);	// Direct Iteration Method
	std::cout << "Solution by DIM:" << std::endl << sol_info.solution;
	std::cout << "Number of iterations: " << sol_info.number_of_it << std::endl << std::endl;

	fout << "Solution by DIM:" << std::endl << sol_info.solution;
	fout << "Number of iterations: " << sol_info.number_of_it << std::endl << std::endl;

	sol_info = Newton_method(start, E, W_of_system, F_of_system);	// Newton's method
	std::cout << "Solution by Newton's method:" << std::endl << sol_info.solution;
	std::cout << "Number of iterations: " << sol_info.number_of_it << std::endl << std::endl;

	fout << "Solution by Newton's method:" << std::endl << sol_info.solution;
	fout << "Number of iterations: " << sol_info.number_of_it << std::endl << std::endl;

	fout.close();
}

//-------------------------------- operator of system -------------------------------

std::vector<double> system_operator(const std::vector<double>& pre_solution) noexcept
{
	std::vector<double> solution(pre_solution.size());
	auto sqr = [](auto elem) -> decltype(elem) {return elem * elem; };
	solution[0] = 0.85 + cos(pre_solution[1]);
	solution[1] = sin(pre_solution[0]) - 1.32;
	return solution;
}

//-------------------------------- Direct Iteration Method --------------------------

info_about_method DIM(std::vector<double> start, double accuracy, std::function<std::vector<double>(const std::vector<double>&)> op) noexcept
{
	auto norm = [](const std::vector<double> v) -> double { return *std::max_element(v.begin(), v.end()); };
	auto pre_sol = start;
	auto sol = op(start);
	info_about_method info;
	info.number_of_it = 1;
	while (norm(sol - pre_sol) >= accuracy)
	{
		pre_sol = sol;
		sol = op(sol);
		info.number_of_it++;
	}
	info.solution = sol;
	return info;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v)
{
	for (const auto& elem : v)
	{
		out << elem << std::endl;
	}
	return out;
}

//-------------------------------- X_n+1 = X_n - W^(-1)(x_n) * F(x_n) ---------------

//-------------------------------- F(x_n) -------------------------------------------

std::vector<double> F_of_system(const std::vector<double>& x) noexcept
{
	std::vector<double> result(x.size());
	result[0] = sin(x[0]) - x[1] - 1.32;
	result[1] = cos(x[1]) - x[0] + 0.85;
	return result;
}

//-------------------------------- W(x_n) -------------------------------------------

mtrx::Square_Matrix<double> W_of_system(const std::vector<double>& x) noexcept
{
	mtrx::Square_Matrix<double> result(x.size());
	result = {
		{ cos(x[0]),	-1			},
		{ -1,			-sin(x[1])	}
	};
	return result;
}

//-------------------------------- Newton's method ----------------------------------

info_about_method Newton_method(std::vector<double> start, double accuracy, std::function<mtrx::Square_Matrix<double>(const std::vector<double>&)> W, std::function<std::vector<double>(const std::vector<double>&)> F) noexcept
{
	auto norm = [](const std::vector<double> v) -> double { return *std::max_element(v.begin(), v.end()); };
	auto pre_sol = start;
	auto sol = start - W(start).inverse() * F(start);
	info_about_method info;
	info.number_of_it = 1;
	while (norm(sol - pre_sol) >= accuracy)
	{
		pre_sol = sol;
		sol = sol - W(sol).inverse() * F(sol);
		info.number_of_it++;
	}
	info.solution = sol;
	return info;
}