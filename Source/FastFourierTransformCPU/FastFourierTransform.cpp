#include "FastFourierTransformCPU/FastFourierTransform.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include "FastFourierTransform.h"

double FFT::magnitude(complex a)
{
	return std::sqrt(a.r * a.r + a.i * a.i);
}

FFT::complex FFT::add(complex a, complex b) {
	return complex(a.r + b.r, a.i + b.i);
}

FFT::complex FFT::add(complex a, double b) {
	return complex(a.r + b, a.i);
}

FFT::complex FFT::add(double a, complex b) {
	return complex(a + b.r, b.i);
}

FFT::complex FFT::mult(complex a, complex b) {
	return complex(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r);
}

FFT::complex FFT::mult(complex a, double b) {
	return complex(a.r * b, a.i * b);
}

FFT::complex FFT::mult(double a, complex b) {
	return complex(b.r * a, b.i * a);
}

FFT::complex FFT::polar(double magnitude, double phase) {
	return FFT::complex(std::cos(phase) * magnitude, std::sin(phase) * magnitude);
}

size_t FFT::reverse_bits(size_t x, int n) {
	size_t result = 0;
	for (int i = 0; i < n; i++, x >>= 1)
		result = (result << 1) | (x & 1U);
	return result;
}

void FFT::fft_radix2(std::vector<FFT::complex>& vec) {
	// Length variables
	size_t n = vec.size();
	int levels = 0;  // Compute levels = floor(log2(n))
	for (size_t temp = n; temp > 1U; temp >>= 1)
		levels++;
	if (static_cast<size_t>(1U) << levels != n)
		throw std::domain_error("Length is not a power of 2");

	std::vector<FFT::complex> expTable(n / 2);
	size_t i;
	for (i = 0; i < n / 2; i++)
		expTable[i] = polar(1.0f, -2 * FFT::PI * i / n);

	for (i = 0; i < n; i++) {
		size_t j = reverse_bits(i, levels);
		if (j > i)
			std::swap(vec[i], vec[j]);
	}

	size_t halfsize;
	size_t tablestep;
	FFT::complex temp;
	size_t size, j, k;
	std::vector<FFT::complex> temp_buffer = vec;

	for (size = 2; size <= n; size *= 2) {
		halfsize = size / 2;
		tablestep = n / size;
		for (i = 0; i < n; i += size) {
			for (j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
				vec[j + halfsize] = add(temp_buffer[j], mult(-1, mult(temp_buffer[j + halfsize], expTable[k])));
				vec[j] = add(temp_buffer[j], mult(temp_buffer[j + halfsize], expTable[k]));
			}
		}
		if (size == n)  // Prevent overflow in 'size *= 2'
			break;
		temp_buffer = vec;
	}
}

void FFT::inverse_fft_radix2(std::vector<FFT::complex>& vec) {

	for (int i = 0; i < vec.size(); i++) {
		//vec[i] = FFT::complex(vec[i].i, vec[i].r);
		vec[i].i *= -1;
	}
	
	FFT::fft_radix2(vec);
	
	for (int i = 0; i < vec.size(); i++) {
		//vec[i] = FFT::complex(vec[i].i, vec[i].r);
		vec[i].i *= -1;
	}

	for (int i = 0; i < vec.size(); i++)
		vec[i] = FFT::complex(vec[i].r / vec.size(), vec[i].i / vec.size());
}

void FFT::fft_shift(std::vector<FFT::complex>& vec)
{
	fft_shift(vec, vec.size() / 2);
}

void FFT::fft_shift(std::vector<FFT::complex>& vec, int shift_amount)
{
	std::vector<FFT::complex> vec_copy = vec;
	for (int i = 0; i < vec.size(); i++) {
		int shifted_index = (i + shift_amount) % vec.size();
		vec[shifted_index] = vec_copy[i];
	}
}

void FFT::inverse_fft_shift(std::vector<FFT::complex>& vec, int shift_amount)
{
	fft_shift(vec, -shift_amount);
}

void FFT::inverse_fft_shift(std::vector<FFT::complex>& vec)
{
	inverse_fft_shift(vec, vec.size() / 2);
}

std::vector<FFT::complex> FFT::read_signal_from_file(const std::string& filepath) {
	std::fstream file(filepath, std::ios::in);
	std::vector<FFT::complex> signal;
	std::string word;
	while (!file.eof()) {
		file >> word;
		signal.push_back(FFT::complex(std::stof(word), 0));
	}
	return signal;
}

std::ostream& operator<<(std::ostream& stream, const FFT::complex& complex) {
	return stream << "(" << complex.r << ", " << complex.i << "i)";
}

