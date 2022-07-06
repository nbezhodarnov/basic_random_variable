#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>

static int y = 0;
static const unsigned long Y_VALUE_LIMIT = (unsigned long)1 << (sizeof(y) * 8 - 1);

// Генератор псевдослучайных чисел (равномерное распределение на отрезке [0, 1])
double rnd() {
	y *= 843314861;
	y += 453816693;
	if (y < 0) {
		y += Y_VALUE_LIMIT;
	}
	return double(y) / double(Y_VALUE_LIMIT - 1);
}

// Аппроксимация квантилей хи-квадрат распределения
double hi_square_distribution_quantile(double alpha, unsigned int degree) {
	if (0 >= alpha || alpha >= 1 || degree == 0) {
		std::cerr << "Wrong alpha or degree value!\n";
		return 0;
	}
	double result = degree, d;
	if (alpha < 0.5) {
		d = -2.0637 * pow(log(1.0 / alpha) - 0.16, 0.4274) + 1.5774;
	} else {
		d = 2.0637 * pow(log(1.0 / (1 - alpha)) - 0.16, 0.4274) - 1.5774;
	}
	result += d * pow(2 * degree, 0.5);
	result += (2.0 / 3.0) * (pow(d, 2) - 1);
	result += d * (pow(d, 2) - 7) / (9.0 * pow(2 * degree, 0.5));
	result += (6 * pow(d, 4) + 14 * pow(d, 2) - 32) / (405.0 * degree);
	result += d * (9 * pow(d, 4) + 256 * pow(d, 2) - 433) / (4860.0 * degree * pow(2 * degree, 0.5));
	return result;
}


int main() {
	const int n = 100;
	std::vector<double> x(n);
	
	// Вычисление числа степеней свободы
	int k = (int)(1 + 3.3 * log10(n));
	std::cout << "Number of degrees of freedom = " << k - 1 << "\n\n";
	
	const double begin = 0, end = 1, probability = 1.0 / (double)k, alpha = 0.9;
	
	// Генерирование случайных чисел
	for (int i = 0; i < n; i++) {
		x[i] = rnd();
	}
	
	// Сортировка по неубыванию
	std::sort(x.begin(), x.end());
	
	// Вывод случайных чисел
	for (int i = 0; i < n; i++) {
		std::cout << i + 1 << ": " << x[i] << '\n';
	}
	std::cout << '\n';
	
	// Вычисление правых границ отрезков, на которые разбивается исходный отрезок (число степеней свободы + 1)
	std::vector<double> interval_ends(k);
	for (int i = 0; i < k; i++) {
		interval_ends[i] = (end - begin) * (i + 1) / (double)k + begin;
	}
	
	// Вычисление количества чисел, которые принадлежат конкретному отрезку
	std::vector<int> observation_numbers(k);
    for (unsigned int i = 0, j = 0; i < n; i++) {
		while (j < interval_ends.size() && x[i] > interval_ends[j]) {
			j++;
		}
		observation_numbers[j]++;
	}
	
	// Вычисление статистического и критического значений
	double statistical_hi = 0, critical_hi = hi_square_distribution_quantile(alpha, k - 1);
	for (int i = 0; i < k; i++) {
		statistical_hi += pow(observation_numbers[i] - n * probability, 2) / (double)(n * probability); 
	}
	std::cout << "Statistical hi^2 = " << statistical_hi << '\n';
	std::cout << "Critical hi^2 (" << alpha << ", " << k - 1 << ") = " << critical_hi << '\n';
	
	// Оценка гипотезы
	if (statistical_hi < critical_hi) {
		std::cout << "There is " << alpha << " probability that the generator obeys a uniform distribution law.\n";
	} else {
		std::cout << "The generator does not obey a uniform distribution law.\n";
	}
	
	return 0;
}
