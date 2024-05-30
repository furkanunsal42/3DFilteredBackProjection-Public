#pragma once

#include <vector>
#include <string>

// fast fourier transform
namespace FFT {

	const double PI = 3.14154965f;

	struct complex {
		complex(double r = 0, double i = 0) :
			r(r), i(i) {}
		double r;
		double i;
	};

	double magnitude(complex a);

	complex add(complex a, complex b);
	complex add(complex a, double b);
	complex add(double a, complex b);

	complex mult(complex a, complex b);
	complex mult(complex a, double b);
	complex mult(double a, complex b);

	complex polar(double magnitude, double phase);

	size_t reverse_bits(size_t x, int n);

	void fft_radix2(std::vector<complex>& vec);
	void inverse_fft_radix2(std::vector<complex >& vec);

	void fft_shift(std::vector<FFT::complex>& vec);
	void fft_shift(std::vector<FFT::complex>& vec, int shift_amount);
	void inverse_fft_shift(std::vector<FFT::complex>& vec, int shift_amount);
	void inverse_fft_shift(std::vector<FFT::complex>& vec);

	std::vector<complex> read_signal_from_file(const std::string& filepath);
}

std::ostream& operator<<(std::ostream& stream, const FFT::complex& complex);